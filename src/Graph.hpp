#ifndef SEMESTRALKA_GRAPH_HPP
#define SEMESTRALKA_GRAPH_HPP

#include <string>
#include <limits>
#include <bitset>
#include <vector>
#include <deque>
#include <omp.h>
#include <mpi.h>
#include "State.hpp"
#include "Solution.hpp"

#define TO_MASTER 0

class Graph {

public:
    explicit Graph(char * argv[]);

    void solve_mpi();

    void solve_data(State stateRecv);

private:
    enum msgTag {
        TASK,
        RESULT,
        DIE
    };

    void loadGraph(const std::string &inputFile);

    void master();

    void sendInitialTasks(std::deque<State> &states, int slaveCount);

    int receiveSolution();

    void sendNextState(std::deque<State> &states, int to);

    void slave();

    std::deque<State> generateStates(State state);

    bool checkState(State & state) const;

    int calculateCutIndividual(std::bitset<Solution::MAX_VERTEX>& setX, int toIndex, int vertex);

    int lowerBound(std::bitset<Solution::MAX_VERTEX>& setX, int index);

    void dfs_par(std::bitset<Solution::MAX_VERTEX> &setX, int currentIndex, int xCount, int cut);

    void dfs_ser(std::bitset<Solution::MAX_VERTEX> &setX, int currentIndex, int xCount, int cut);

    void printSet(std::ostream &s, bool setX) const;

    void master_work(State & stateRecv);

    void printEdges(std::ostream &s) const;

    State addToSet(State state, bool toX);

    friend std::ostream& operator<<(std::ostream& s, const Graph& g);

    constexpr static const int m_stateSize = sizeof(class State);
    constexpr static const int m_solutionSize = sizeof(class Solution);
    int m_adj[Solution::MAX_VERTEX][Solution::MAX_VERTEX];
    Solution m_solution;
    int m_vertexCount;
    int m_a;
    int m_maxBFSDepth;
    size_t m_stateCount;
    int m_maxDFSParDepth;
    int m_procCount;
    int m_rank;
    int m_threadCount;
    std::string m_filename;
};


#endif //SEMESTRALKA_GRAPH_HPP
