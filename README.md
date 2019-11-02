# kt_triangle

Code for investigating properties of triangle diagrams using dispersive techniques in the context of Khuri-Treiman.

Requires ROOT (tested with version 6.17).

### Executables
To build any example executable, for example `test.cpp`, use:
```
mkdir build
cd build
cmake ..
make test
````

The `plot_triangle.cpp` executable simply evaluates the triangle functions numerically for a given set of masses as a function of s in the Feynman representation and dispersion relation representation. These should match. Produces .dat files containing s and real, imaginary, and absolute values of the triangle function calculated with each method.
