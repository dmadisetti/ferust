#!/usr/bin/python3
""" generate.py
This script computes key values in finite element, including the stiffness
matrix, the mass matrix, the projection vector, and interplation functions
symbolically, and then compiles the result down into an optimized matlab
form. Note, the performance gains aren't as great as I hoped, but the idea was
fun.
"""
import symengine as sympy
from symengine import var, expand, cse
import numpy as np
import pystache
from itertools import product
import re
import time
import sympy as sympy2
from multiprocessing import Pool
import pdb

# Define loose variables
(  # General system variables
    n, e, h, f, l, u, h,
    # Loose variables for X_0 -> X_9
    A, B, C, D, E, F, G, H, I,
    # Loose Variables for Y_0 -> Y_9
    R, S, T, U, V, W, X, Y, Z,
) = var("n, e, h, f, l, u, h,"
        "A, B, C, D, E, F, G, H, I,"
        "R, S, T, U, V, W, X, Y, Z")


# Shape functions
N = lambda E, N: ((1 + e * E) * (1 + n * N)) / 2
def quad():
    N1 = N(-1, -1) / 2
    N2 = N(1, -1) / 2
    N3 = N(1, 1) / 2
    N4 = N(-1, 1) / 2
    return [N1, N2, N3, N4]


def bubble(N5=True, N6=True, N7=True, N8=True, N9=True):
    (N1, N2, N3, N4) = quad()
    N9 = 2 * int(N9) * N(-e, -n)

    N1 -= N9 / 4
    N2 -= N9 / 4
    N3 -= N9 / 4
    N4 -= N9 / 4

    N5 = int(N5) * (N(-e, -1) - N9 / 2)
    N6 = int(N6) * (N(1, -n) - N9 / 2)
    N7 = int(N7) * (N(-e, 1) - N9 / 2)
    N8 = int(N8) * (N(-1, -n) - N9 / 2)

    N1 -= (N5 + N8) / 2
    N2 -= (N5 + N6) / 2
    N3 -= (N6 + N7) / 2
    N4 -= (N7 + N8) / 2
    return [N1, N2, N3, N4, N5, N6, N7, N8, N9]


def serindipity(N5=True, N6=True, N7=True, N8=True):
    Ns = bubble(N5=N5, N6=N6, N7=N7, N8=N8, N9=False)
    return Ns[:8]

# Global definitions
combination = lambda *a: list(product(a, a))
shapes = {4: quad, 8: serindipity, 9: bubble}
integration_points = {
    4:
    combination(-1 / sympy.sqrt(3), 1 / sympy.sqrt(3)),
    9:
    combination(-sympy.sqrt(3) / sympy.sqrt(5),
                sympy.sqrt(3) / sympy.sqrt(5), 0)
}
integration_weights = {
    4: combination(1, 1),
    9: combination(sympy.S(5) / 9,
                   sympy.S(5) / 9,
                   sympy.S(8) / 9)
}
weight_lookup = {4: 4, 8: 9, 9: 9}

x0 = [A, B, C, D, E, F, G, H, I]
y0 = [R, S, T, U, V, W, X, Y, Z]
disp = [tuple(var(f"d{i}{['X', 'Y'][j]}") for j in range(2)) for i in range(9)]
placeholder = [var(f"d{i}") for i in range(9)]

# We only do the plane strain case, since a different l will be passed down.
D = [[l + 2 * u, l, 0], [l, l + 2 * u, 0], [0, 0, u]]


def B(N, x, y, J):
    Nx = N.diff(e) * y.diff(n) - N.diff(n) * y.diff(e)
    Ny = -N.diff(e) * x.diff(n) + N.diff(n) * x.diff(e)
    return np.array([[Nx, 0, Ny], [0, Ny, Nx]]).T / J


def format_code(code,
                convert=True,
                trim=2,
                padding=13,
                cols=50,
                initial_offset=4,
                sep=" ...\n    "):
    """Massages generated code into a usable form. Lots of strings hacks, but it
    seems to work!"""
    inner_width = cols - padding
    initial_width = inner_width - initial_offset

    if convert:
        code = sympy2.octave_code(code)[trim:-trim]

    K = code
    code = code.replace("{", "[").replace("}", "]").replace(".*", "*").replace(
        "], [", "; ").replace("**", "^").replace(".^", "^")
    #code = re.sub("\s", "", code)
    code = re.sub("sqrt\(", "&", code)
    MATCH_REGEX = "([dx]?[0-9e\.-]+[XY]?)"
    matches = re.findall(MATCH_REGEX, code)
    code = re.sub(MATCH_REGEX, "?", code)

    code = code[:initial_width] + sep + sep.join([
        code[i:i + inner_width]
        for i in range(initial_width, len(code), inner_width)
    ])
    return code.replace('?', "{}").format(*matches).replace("&", "sqrt(")

def horner(K):
    # Something I was prototyping but eventually scrapped.
    # symbols = list(sorted(K.free_symbols))
    # k = expand(expand(K))
    # P = np.array([[
    #     int(list(ex.subs(
    #         symb, 2).as_coefficients_dict().values())[0]).bit_length() - 1
    #     for symb in symbols
    # ] for ex in k.as_coefficients_dict()])
    return K


def factor(k):
    # Ideally reduces the equation, but this is very slow.
    if isinstance(k, int):
        return k
    return k
    #frac = k.as_numer_denom()
    #frac = list(map(horner, frac))
    #return frac[0] / frac[1]

def symbolic_interpolation(nodes):
    """Exports our interpolation functions."""
    Ns = shapes[nodes]()

    x = sum([n * x_ for n, x_ in zip(Ns, x0[:nodes])])
    y = sum([n * y_ for n, y_ in zip(Ns, y0[:nodes])])

    s = sum([N * interp for N, interp in zip(Ns, placeholder[:nodes])])
    s = factor(s)
    return """
    fn = @(n, e) {};""".format(
        format_code(sympy2.octave_code(s), convert=False, cols=90))


def symbolic_equations(nodes, weights=0, dim=2):
    """Computes values for local_values and local_projection."""
    stress_components = {1: 1, 2: 3, 3: 6}[dim]
    if not weights:
        weights = weight_lookup[nodes]

    Ns = shapes[nodes]()

    x = sum([n * x_ for n, x_ in zip(Ns, x0[:nodes])])
    y = sum([n * y_ for n, y_ in zip(Ns, y0[:nodes])])

    J = (sympy.diff(x, e) * sympy.diff(y, n) -
         sympy.diff(x, n) * sympy.diff(y, e))

    stress = np.reshape(
        sum([
            np.matmul(np.matmul(D, B(N, x, y, J)), (dx, dy))
            for N, (dx, dy) in zip(Ns, disp[:nodes])
        ]), (stress_components, 1))

    # a list of Our integration points and weightings.
    X = integration_points[weights]
    W = integration_weights[weights]

    # Init K and F to 0
    K = np.zeros((nodes * dim, nodes * dim), dtype="object")
    F = np.zeros((nodes * dim, 1), dtype="object")
    # Init value for our stress calculations
    M = np.zeros((nodes * stress_components, nodes * stress_components),
                 dtype="object")
    P = np.zeros((nodes * stress_components, 1), dtype="object")

    subs = lambda Xe, Xn: lambda k: list(
        map(lambda k: k.subs({
            e: Xe,
            n: Xn
        }), k))
    for i, Na in enumerate(Ns):
        for j, Nb in enumerate(Ns):
            # symmetric
            if i > j:
                continue

            print(i, j)

            k = np.matmul(np.matmul(B(Na, x, y, J).T, D), B(Nb, x, y, J)) * J
            k = np.sum([
                Wx * Wy * np.array(list(map(subs(Xe, Xn), k)))
                for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
            ],
                       axis=0)
            K[dim * i:dim * i + dim, dim * j:dim * j + dim] += k

            # So we end up with just URH for any dimension
            if i == j:
                for clear_i in range(dim):
                    for clear_j in range(dim):
                        if clear_i > clear_j:
                            K[dim * i + clear_i, dim * i + clear_j] = 0

            # Calculate Mass Matrix contributions
            m = Na * Nb * J
            m = np.sum([
                Wx * Wy * m.subs({
                    e: Xe,
                    n: Xn
                }) for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
            ],
                       axis=0) * np.eye(stress_components)
            M[stress_components * i:stress_components * i +
              stress_components, stress_components * j:stress_components * j +
              stress_components] += m
            # So we end up with just URH for any dimension
            if i == j:
                for clear_i in range(stress_components):
                    for clear_j in range(stress_components):
                        if clear_i > clear_j:
                            M[stress_components * i +
                              clear_i, stress_components * i + clear_j] = 0

        # We can concurrently solve for the base F component.
        # We only set the body force from the y component
        temp = Na * f * J
        F[1 + i * dim] = sum([
            Wx * Wy * temp.subs({
                e: Xe,
                n: Xn
            }) for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
        ])

        # We can also set up our projection vector
        p = Na * stress * J
        P[stress_components * i:stress_components * i +
          stress_components] = sum([
              Wx * Wy * np.array(list(map(subs(Xe, Xn), p)))
              for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
          ])

    # We put it in horner form since it is supposed to be one of the best forms
    # for polynomial evaluation.
    K = list(map(factor, np.reshape(K, (1, -1))[0]))
    M = list(map(factor, np.reshape(M, (1, -1))[0]))
    F = np.array(
        list(map(lambda f: list(map(sympy2.polys.polyfuncs.horner, f)),
                 F))).T[0]
    P = list(map(factor, np.reshape(P, (1, -1))[0]))
    precompute, reduced = cse(np.concatenate((K, M, F)))
    precomputeP, reducedP = cse(P)

    # Reshape our square matrices
    F = reduced[len(K) + len(M):]
    M = np.reshape(reduced[len(K):-len(F)],
                   (nodes * stress_components, nodes * stress_components))
    K = np.reshape(reduced[:len(K)], (nodes * dim, nodes * dim))

    # Generate the CSE components.
    code = ";".join([
        sympy2.octave_code(v) + "=" + sympy2.octave_code(k)
        for (v, k) in precompute
    ])
    codeP = ";".join([
        sympy2.octave_code(v) + "=" + sympy2.octave_code(k)
        for (v, k) in precomputeP
    ])
    return """
    {};
    K=[{}];
    F=[{}]\';
    M=[{}];""".format(format_code(code, convert=False), format_code(K),
                      format_code(F, trim=1), format_code(M)), """
    {};
    P=[{}];""".format(format_code(codeP, convert=False),
                      format_code(reducedP, trim=1))


def main():
    """Generate all the files needed."""
    renderer = pystache.Renderer()

    quad = symbolic_interpolation(4)
    serindipity = symbolic_interpolation(8)
    bubble = symbolic_interpolation(9)
    with open(f"matlab/local_interpolation.m", 'w') as file:
        print(renderer.render_path("templates/local_interpolation.m.tmpl", {
            "quad": quad,
            "serindipity": serindipity,
            "bubble": bubble,
        }),
              file=file)
    for weights in [4, 9]:
        quad, qp = symbolic_equations(4, weights=weights)
        serindipity, qs = symbolic_equations(8, weights=weights)
        bubble, qb = symbolic_equations(9, weights=weights)
        # Generate our stiffness equations.
        with open(f"matlab/local_values_{weights}.m", 'w') as file:
            print(renderer.render_path(
                "templates/local_values.m.tmpl", {
                    "integration": weights,
                    "quad": quad,
                    "serindipity": serindipity,
                    "bubble": bubble,
                }),
                  file=file)
        with open(f"matlab/local_projection_{weights}.m", 'w') as file:
            print(renderer.render_path(
                "templates/local_projection.m.tmpl", {
                    "integration": weights,
                    "quad": qp,
                    "serindipity": qs,
                    "bubble": qb,
                }),
                  file=file)


if __name__ == "__main__":
    main()
