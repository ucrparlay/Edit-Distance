## Edit-Distance

Efficient parallel output-sensitive edit-distance implementation.

Authors: ucrparlay

### Requirements

CMake >= 3.15

### Usage

Clone the repository with submodules
```
$ git clone --recurse-submodules git@github.com:ucrparlay/Edit-Distance.git
$ cd Edit-Distance
```

Build & Run
```
$ mkdir build && cd build
$ cmake ..
$ make
$ ./test_framework -1 1000 10
```

You can specify the algorithm by the first argument of test_framework. "-1" means run all algorithms.
