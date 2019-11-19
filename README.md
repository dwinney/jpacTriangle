# kt_triangle

Code for investigating properties of triangle diagrams using dispersive techniques in the context of Khuri-Treiman.

Requires ROOT (tested with version 6.17).

The primary class of interest is the `triangle` class which calculates the Feynman rescattering diagram as a function of center-of-mass energy for a given set of fixed masses (two internal propagating masses and two external masses) and a given parameterization of the left-hand cut amplitude.

The LHC amplitudes can be user-implemented by inheriting from the abstract `lefthand_cut` class. All that is needed is a way to evaluate the discontinuity across the cut in a dispersive representation. For example:
```
 class my_Model : public lefthand_cut
 {
   my_Model( my_inputs )
   {};

  {
    // All other code...
  }

   complex<double> eval(double s)
   {
     // code to evaluate the whole propagator
   };

   complex<double> disc(double s)
   {
     // code to evaluate the discontinuity
   };
 };
```

Then you can define a `triangle` amplitude to evaluate as:
```
my_Model LHC(my_Inputs);

triangle my_triangle(&LHC);

// Set masses
triangle.set_externalMasses(m1, m2);
triangle set_internalMasses(m3, m4);

// Evaluate
triangle.eval_feynman(s);
```

### Executables
To build any example codes in the /executable/ directly, for example `test.cpp`, use:
```
mkdir build
cd build
cmake ..
make test
````
