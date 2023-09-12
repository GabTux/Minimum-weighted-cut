#include <iostream>
#include <chrono>
#include "Graph.hpp"


int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cout << "Invalid arguments" << std::endl;
        return EXIT_FAILURE;
    }

    static_assert(std::is_trivially_copyable<State>::value, "State must be trivially-copyable");

    int rank;
    int provided, required = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided < required) {
        std::cout << "Bad MPI model" << std::endl;
        return EXIT_FAILURE;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        auto start = std::chrono::steady_clock::now();
        Graph graph(argv);
        graph.solve_mpi();
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration = end-start;
        std::cout << "TIME: " << duration.count() << "s" << std::endl;
        std::cout << graph << std::endl;
    } else {
        Graph graph(argv);
        graph.solve_mpi();
    }

    MPI_Finalize();
    return 0;
}
