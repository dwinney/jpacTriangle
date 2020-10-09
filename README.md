# kt_triangle

Code for investigating properties of triangle diagrams using dispersive techniques in the context of Khuri-Treiman (see for example [[1]](https://arxiv.org/abs/1910.03107)).

<p align="center">
  <img width="500" src="./doc/tri_diagram.png">
</p>
The goal is the evaluation of the diagram above with arbitrary external spin-<em>J</em> of the decaying particle, spin-projection-<em>j</em> in the s-channel and exchanged particle of spin-<em>j'</em> in the t-channel (as well as arbitrary helicity projections <em>λ</em> and <em>λ'</em>).

## USAGE

Requires [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed. And [BOOST](https://www.boost.org/) C++ libraries. The later is included for adaptive integration of a complex-valued function (see [here](https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/gauss_kronrod.html))

This repo uses `git submodules` so when cloning make sure to use:
```
git clone --rescursive https://github.com/dwinney/jpacTriangle.git
```

To build linkable libraries `jpacTriangle` and `jpacStyle` use:
```bash
mkdir build && cd build
cmake ..
cmake --build . --target install
````

## REFERENCES
* [1] "Khuri-Treiman equations for 3π decays of particles with spin" JPAC Collaboration [[arXiv:1910.03107]](https://arxiv.org/abs/1910.03107)
