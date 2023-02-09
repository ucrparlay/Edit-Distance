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
$ ./test_framework -1 1000 10    # For randomly generated sequences.
$ ./test_framework_real -i <algorithm id> -f1 <path to file1> -f2 <path to file2>    # For real world data
```

You can specify the algorithm by the first argument of test_framework. "-1" means run all algorithms. You can also specify the algorithm by selecting <algorithm id>, and any two text file paths <path to file1>, <path to file2>.

