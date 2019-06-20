# oriented-particles
Oriented Particle System simulations

[![MainWindow](https://raw.githubusercontent.com/amit112amit/oriented-particles/master/src/Gui/YouTubeLink.png)](https://youtu.be/IOLdgTWHIAM)

Description:
1. Creates a wrapper to l-BFGS-b Fortran library.
2. Implements Oriented Particle System (OPS) energy and jacobian for minimization.
3. Implements Brownian dynamics on a curved surface using OPS.

Dependencies:
1. VTK
2. Eigen
3. lbfgsb
4. SHTOOLS

## Building on Windows 10 with MSYS2 or MinGW-W64:
1. We must manually specify the path for `qwt` headers and the library. Look at the end of `CMakeLists.txt` in the root directory.
2. We must manually specify the path for `shtools` module and library. e.g. `SHTOOLS_LIBRARY=C:/msys64/usr/local/lib/libSHTOOLS.a` and `SHTOOLS_INCLUDE_DIR=C:/msys64/usr/local/include`. SHTOOLS project does not use CMake. It only provides a Makefile. Hence, it is easiest to build it directly under MSYS2.

References:
1. Szeliski, Richard and Tonnesen, D. (1992). Surface Modeling with Oriented Particle Systems. Siggraph ’92, 26(2), 160. https://doi.org/10.1017/CBO9781107415324.004
2. Liu, D.C. & Nocedal, J. Mathematical Programming (1989) 45: 503. https://doi.org/10.1007/BF01589116
3. SHTOOLS https://shtools.oca.eu/shtools/
4. https://github.com/wsklug/voom
