#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include "Graph.hpp"
#include "State.hpp"


Graph::Graph(char **argv) :
m_solution(),
m_filename(argv[1])
{
    if (!std::filesystem::exists(m_filename))
        throw std::runtime_error("WRONG FILE");
    std::istringstream ss(argv[2]);
    ss >> m_a;
    ss = std::istringstream(argv[3]);
    ss >> m_threadCount;
    omp_set_num_threads(m_threadCount);
    loadGraph(m_filename);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_procCount);
    m_maxBFSDepth = static_cast<int>(m_vertexCount * 0.7);
    m_maxDFSParDepth = static_cast<int>(m_vertexCount * 0.6);
    m_stateCount = (m_procCount-1)*2;
}


void Graph::loadGraph(const std::string& inputFile) {
    std::ifstream inputStream(inputFile);
    std::string line;
    std::getline(inputStream, line);
    std::istringstream issFirstLine(line);
    issFirstLine >> m_vertexCount;
    for (int rowAdj = 0; rowAdj < m_vertexCount; ++rowAdj) {
        std::getline(inputStream, line);
        std::istringstream iss(line);
        for (int columnAdj = 0; columnAdj < m_vertexCount; ++columnAdj)
            iss >> m_adj[rowAdj][columnAdj];
    }
}

int Graph::calculateCutIndividual(std::bitset<Solution::MAX_VERTEX>& setX, int toIndex, int vertex) {
    // assumes that vertices < toIndex have been placed to sets
    // calculate cost cut only for vertex
    int cut = 0;
    for (int n = 0; n < toIndex; n++) {
        if (m_adj[vertex][n]) {
            if (setX[vertex] != setX[n]) {
                cut += m_adj[vertex][n];
            }
        }
    }
    return cut;
}

int Graph::lowerBound(std::bitset<Solution::MAX_VERTEX>& setX, int index) {
    if (m_solution.m_minCut == std::numeric_limits<int>::max())
        return 0;

    int lowerBound = 0;
    for (int i = index; i < m_vertexCount; ++i) {
        setX[i] = true;
        int x = calculateCutIndividual(setX, index, i);
        setX[i] = false;
        int y = calculateCutIndividual(setX, index, i);
        lowerBound += std::min(x, y);
    }

    return lowerBound;
}

void Graph::printSet(std::ostream &s, bool setX) const {
    for (int i = 0; i < m_vertexCount; ++i) {
        if (m_solution.m_result[i] && setX) {
            s << i << " ";
        } else if (!m_solution.m_result[i] && !setX) {
            s << i << " ";
        }
    }
}

void Graph::printEdges(std::ostream &s) const {
    for (int a = 0; a < m_vertexCount; ++a) {
        for (int b = a; b < m_vertexCount; ++b) {
            if (m_adj[a][b]) { // neighbours?
                if (m_solution.m_result[a] != m_solution.m_result[b]) { // different sets?
                    s << a << "--" << m_adj[a][b] << "--" << b << "  ";
                }
            }
        }
    }
}

std::ostream &operator<<(std::ostream &s, const Graph &g) {
    s << g.m_filename << " " << g.m_a << std::endl;
    s << "Threads: " << g.m_threadCount << std::endl;
    s << "MinCut: " << g.m_solution.m_minCut << std::endl;

    s << "X: ";
    g.printSet(s, true);

    s << std::endl << "Y: ";
    g.printSet(s, false);

    s << std::endl << "Edges: ";
    g.printEdges(s);

    s << std::endl;
    return s;
}

void Graph::solve_mpi() {
    if (m_rank == 0)
        master();
    else
        slave();
}

void Graph::dfs_par(std::bitset<Solution::MAX_VERTEX> &setX, int currentIndex, int xCount, int cut) {
    // invalid solutions
    if (((m_a - xCount) > (m_vertexCount - currentIndex)) // not enough vertices left for X
        || (cut >= m_solution.m_minCut)
        || ((cut+ lowerBound(setX, currentIndex)) >= m_solution.m_minCut)
            )
    {
        return;
    }

    // X is full, assign rest to the Y
    if (xCount == m_a) {
        while (currentIndex < m_vertexCount) {
            setX[currentIndex] = false;
            cut += calculateCutIndividual(setX, currentIndex, currentIndex);
            currentIndex++;
        }
    }

    if (currentIndex == m_vertexCount) {
        // All vertices have been assigned, and we have valid solution (checks above)
        if (cut < m_solution.m_minCut) {
            #pragma omp critical
            if (cut < m_solution.m_minCut) {
                m_solution = Solution(cut, setX);
            }
        }
        return;
    }

    // assign vertex to X
    setX[currentIndex] = true;
    if (currentIndex < m_maxDFSParDepth) {
        #pragma omp task default(none) firstprivate(setX) \
        shared(currentIndex, xCount, cut, m_solution)
        dfs_par(setX, currentIndex + 1, xCount + 1, cut + calculateCutIndividual(setX, currentIndex, currentIndex));
    } else {
        dfs_ser(setX, currentIndex + 1, xCount + 1, cut + calculateCutIndividual(setX, currentIndex, currentIndex));
    }

    // assign vertex to Y
    setX[currentIndex] = false;
    if (currentIndex < m_maxDFSParDepth) {
        #pragma omp task default(none) firstprivate(setX) \
        shared(currentIndex, xCount, cut, m_solution)
        dfs_par(setX, currentIndex + 1, xCount, cut + calculateCutIndividual(setX, currentIndex, currentIndex));
    } else {
        dfs_ser(setX, currentIndex + 1, xCount, cut + calculateCutIndividual(setX, currentIndex, currentIndex));
    }
}


void Graph::dfs_ser(std::bitset<Solution::MAX_VERTEX> &setX, int currentIndex, int xCount, int cut) {
    // invalid solutions
    if ((xCount > m_a)
        || ((m_a - xCount) > (m_vertexCount - currentIndex)) // not enough vertices left for X
        || (cut >= m_solution.m_minCut)
        || ((cut+ lowerBound(setX, currentIndex)) >= m_solution.m_minCut)
        )
    {
        return;
    }

    // X is full, assign rest to the Y
    if (xCount == m_a) {
        while (currentIndex < m_vertexCount) {
            setX[currentIndex] = false;
            cut += calculateCutIndividual(setX, currentIndex, currentIndex);
            currentIndex++;
        }
    }

    if (currentIndex == m_vertexCount) {
        // All vertices have been assigned, and we have valid solution (checks above)
        if (cut < m_solution.m_minCut) {
            #pragma omp critical
            if (cut < m_solution.m_minCut) {
                m_solution = Solution(cut, setX);
            }
        }
        return;
    }

    // assign vertex to X
    setX[currentIndex] = true;
    dfs_ser(setX, currentIndex + 1, xCount + 1, cut + calculateCutIndividual(setX, currentIndex, currentIndex));

    // assign vertex to Y
    setX[currentIndex] = false;
    dfs_ser(setX, currentIndex + 1, xCount, cut + calculateCutIndividual(setX, currentIndex, currentIndex));
}


std::deque<State> Graph::generateStates(State state) {
    std::deque<State> qStates;
    qStates.emplace_front(state);

    while (qStates.size() < m_stateCount) {
        State currentState = qStates.back();
        if (currentState.currentIndex >= m_maxBFSDepth)
            break;
        qStates.pop_back();
        State xState = addToSet(currentState, true);
        State yState = addToSet(currentState, false);
        if (checkState(xState))
            qStates.emplace_front(xState);
        if (checkState(yState))
            qStates.emplace_front(yState);
    }
    return qStates;
}

State Graph::addToSet(State state, bool toX) {
    state.setX[state.currentIndex] = toX;
    state.cut += calculateCutIndividual(state.setX, state.currentIndex, state.currentIndex);
    state.currentIndex++;
    if (toX)
        state.xCount++;
    return state;
}

bool Graph::checkState(State & state) const {
    if ((state.xCount > m_a)
        || ((m_a - state.xCount) > (m_vertexCount - state.currentIndex))
        )
    {
        return false;
    }
    return true;
}

void Graph::master() {
    // generate states for slaves
    std::deque<State> states = generateStates(State());
    int slaveCount = m_procCount - 1;

    sendInitialTasks(states, slaveCount);

    if (!states.empty()) {
        State state = states.back();
        states.pop_back();
        master_work(state);
    }

    while (!states.empty()) {
        int from = receiveSolution();
        sendNextState(states, from);
    }

    // collect results and kill slaves
    while (slaveCount > 0) {
        char dummy = 0x00;
        int from = receiveSolution();
        MPI_Send(&dummy, 1, MPI_BYTE, from, msgTag::DIE, MPI_COMM_WORLD);
        slaveCount--;
    }
}

void Graph::slave() {
    MPI_Status status;
    State stateRecv;

    while (true) {
        MPI_Recv(&stateRecv, m_stateSize, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == msgTag::DIE) {
            break;
        } else {
            if (stateRecv.currentBest < m_solution) {
                m_solution = stateRecv.currentBest;
            }
            //solve_data(State());
            #pragma omp parallel
            {
                #pragma omp single
                dfs_par(stateRecv.setX, stateRecv.currentIndex, stateRecv.xCount, stateRecv.cut);
            }
            MPI_Send(&m_solution, m_solutionSize, MPI_PACKED, TO_MASTER, msgTag::RESULT, MPI_COMM_WORLD);
        }
    }
}

void Graph::sendInitialTasks(std::deque<State> &states, int slaveCount) {
    for (int i = 1; i <= slaveCount && !states.empty(); ++i) {
        sendNextState(states, i);
    }
}

int Graph::receiveSolution() {
    MPI_Status status;
    Solution solutionRecv;
    MPI_Recv(&solutionRecv, m_solutionSize, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (solutionRecv < m_solution) {
        m_solution = solutionRecv;
    }
    return status.MPI_SOURCE;
}

void Graph::sendNextState(std::deque<State> &states, int to) {
    State state = states.back();
    states.pop_back();
    state.currentBest = m_solution;
    MPI_Send(&state, m_stateSize, MPI_PACKED, to, msgTag::TASK, MPI_COMM_WORLD);
}

void Graph::solve_data(State stateRecv) {
    std::deque<State> qStates(generateStates(stateRecv));
    std::vector<State> states(qStates.begin(), qStates.end());
    //std::sort(states.begin(), states.end());

    #pragma omp parallel for default(shared) schedule(dynamic)
    for (auto & state : states) {
        dfs_par(state.setX, state.currentIndex, state.xCount, state.cut);
    }
}

void Graph::master_work(State &stateRecv) {
    #pragma omp parallel
    {
        #pragma omp single
        dfs_par(stateRecv.setX, stateRecv.currentIndex, stateRecv.xCount, stateRecv.cut);
    }
}
