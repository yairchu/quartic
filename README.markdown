# Quartic

Solve Quartic equations (4th degree polynomials) in C, Haskell, and Python.

License: BSD.

(Also solves cubic, quadratic, and linear equations)

The implementations are very close in their structure/design.
Naturally, the Haskell and Python implementations are much shorter than the C implementation.

## C implementation

There's one implementation for complex numbers, and another for real numbers (`double`), which only finds real roots.

## Haskell implementation

This implementation is generic and should work for finite fields too (I think).

### Note about Cubic equations

The solving of cubic equations only finds one root.
Finding all three roots for `Complex Double`s isn't hard but I wanted to keep the code generic.

For that a type-class providing primitive roots of unity would be required.
The class providing the primitive square root of unity, `(-1)`, exists in Haskell and it is `Floating`, but here we would need a primitive cubic root of unity..

## Python implementation

Numpy already has a built-in roots function, so this implementation isn't necessary.
It was only made for comparison with the others.
Also consider https://github.com/NKrvavica/fqs if you require a faster implementation than numpy.