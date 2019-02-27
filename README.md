# topopt99
A simple, easily hackable topopt code in C/C++

This code is directly inspired from the fantastic introductory paper:

*A 99 line topology optimization code written in MATLAB*<br>
Structural and Multidisciplinary Optimization 21(2), 2001, pp. 120-127<br>
by Ole Sigmund.

http://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB

A long time ago I re-implemented it in C/C++. I thought it would make sense to share it for educational purposes, as I found it very helpful as a starting point to play with topology optimization and build an intuition on how it actually works. This updated version uses LibSL-small and Eigen.<br>

And yes, it takes more than 99 lines of C/C++ ;-)

### How to
- Clone repo
- Clone submodules (git submodule init; git submodule update)
- Compile with CMake
- Run, tga images of the structure are output in the current directory at each iteration

