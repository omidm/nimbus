# Evaluation and Profiling


## How to profile code using `perf`

First, compile the code using `-fno-omit-frame-pointer -ggdb -pg` flags.
Uncomment the line in `Makeinclude` file at the Nimbus root that is for
profiling with `perf`. Then clean the previous build and build again:

    $ make clean-hard
    $ make

Then, run the executables with `perf`:

    $ perf record -g <executable> <args>

It dumps `per.data` file, that can be processed. To get the call stacks based
on the caller trace.

    $ perf report -g fractal,0.5,caller -i <dumped-file>

To produce `FlameGraph`, download the repo from `github`. Run the following
command in the directory that has `perf.data`:

    $ perf script | ./stackcollapse-perf.pl | ./flamegraph.pl > result.svg

##  How to remove DBG in compile time

1. Uncomment `-D_NIMBUS_NO_DBG` line from `Makeinclude` file at the root.

2. Uncomment `add_definitions(-D_NIMBUS_NO_DBG)`
   from `application/water_multiple/CMakeLists.txt`

3. Comment out `ADD_DEFINITIONS (-Werror)`
   from `application/physbam-lib/Scripts/CMake/Compiler_Flags.cmake`

4. Make sure that you `make clean-hard` and `make` to rebuild.



