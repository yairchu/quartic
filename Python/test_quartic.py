# Run with: python -m pytest

from hypothesis import given, strategies as st

import quartic

reasonable_complex = st.complex_numbers(
    max_magnitude=100, min_magnitude=0.01, allow_nan=False
)
complex_2_to_4 = st.lists(reasonable_complex, min_size=2, max_size=4)


@given(complex_2_to_4)
def test_solve_poly(poly):
    if poly[0] == 0 or poly[-1] == 0:
        # Skip testing too-small polynomials
        return
    roots = quartic.stable_solve_poly(poly)
    if not roots:
        for x in poly:
            assert x == 0
        return
    for root in roots:
        x = evaluate_poly(poly, root)
        assert abs(x) < 0.01, f"Result at root {root} too large: {x}"


@given(complex_2_to_4)
def test_from_roots(vals):
    poly = poly_from_roots(vals)
    roots = quartic.solve_poly(poly)
    elem_check(vals, roots)
    if len(vals) != 3:
        elem_check(roots, vals)


def elem_check(expected, got):
    "Verify the elements from given set appear up to a small error within another set"
    for x in got:
        for y in expected:
            if abs(x - y) < 0.01:
                break
        else:
            assert False, f"Value {x} not found in {expected}"


def poly_from_roots(roots):
    result = [1]
    for root in roots:
        result = convolve(result, [root, -1])
    return result


def convolve(xs, ys):
    result = [0] * (len(xs) + len(ys) - 1)
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            result[ix + iy] += x * y
    return result


# TODO: Also implement propsRealQuadratic test


def evaluate_poly(poly, x):
    result = 0
    p = 1
    for coef in poly:
        result += p * coef
        p *= x
    return result
