# Minimum-weighted-cut
This C++ project leverages MPI (Message Passing Interface) and OpenMP to address the challenge of finding the minimum edge cut of a weighted graph. MPI and OpenMP are employed to harness parallelism and enable potential execution on cluster environments.

### Usage
1) clone:
```
git clone https://github.com/GabTux/Minimum-weighted-cut
```
2) build:
```
cd Minimum-weighted-cut
cmake CMakeLists.txt
make
```
3) run:
```
mpiexec -n <num of processes> ./minCut <path_to_input.txt> <cut size> <num of threads for one process>
```

Example:
```
mpiexec -n 4 ./minCut ../input/graf_40_8.txt 20 14
```
