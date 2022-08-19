# MaCh3_DUNE

## Dependencies

- MaCh3 Core
- C++17
- CMake 3.17+


## Build

### Quick-start

```
$ git clone mach3-software/MaCh3
$ cd MaCh3; git checkout feature/removelibconfig_addyamlcpp;
$ cd ../;
$ git clone luketpickering/MaCh3_DUNE
$ cd MaCh3_DUNE; git checkout feature/CPMCore;
$ mkdir build; cd build;
$ cmake .. -DCPM_MaCh3_SOURCE=../../MaCh3
$ make -j 10
```

### MaCh3 Core

This build uses CPM to include configure-, build- and run-time dependencies. CPM will attempt to pull down MaCh3 core for you, but if MaCh3 core is in a private repository, it will require you to manually put in your github credentials. Instead you should clone `mach3-software/MaCh3` yourself and point CPM to the repository like by passing the `-DCPM_MaCh3_SOURCE=/path/to/MaCh3` flag at configure time. **N.B.** You do not need to build MaCh3 manually, MaCh3_DUNE will include it as a subproject and build it within the MaCh3_DUNE build tree. 

### Optional Falgs

```
$ cmake ../MaCh3_DUNE -DCPU_ONLY=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples 
#CUDA_SAMPLES not necessary if using CPU_ONLY=ON
$ make
```