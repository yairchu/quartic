# Quartic

Solve Quartic equations (4th degree polynomials) in C and Haskell.

License: BSD.

(Also solves cubic, quadratic, and linear equations)

The two implementations are very close in their structure/design.
Naturally, the Haskell implementation is shorter.

## C implementation

There's one implementation for complex numbers, and another for real numbers (`double`), which only finds real roots.

## Haskell implementation

This implementation is generic and should work for finite fields too (I think).

### Note about Cubic equations

The solving of cubic equations only finds one root.
Finding all three roots for `Complex Double`s isn't hard but I wanted to keep the code generic.

For that a type-class providing primitive roots of unity would be required.
The class providing the primitive square root of unity, `(-1)`, exists in Haskell and it is `Floating`, but here we would need a primitive cubic root of unity..
