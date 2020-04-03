#!/usr/bin/python3
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
    # Loose variables for X_0 -> X_8
    A, B, C, D, E, F, G, H,
    # Loose Variables for Y_0 -> Y_8
    S, T, U, V, W, X, Y, Z
) = var("n, e, h, f, l, u, h,"
             "A, B, C, D, E, F, G, H, S,"
             "T, U, V, W, X, Y, Z")

# using the shape function convention in #5
N = lambda E, N: ((1 + e * E) * (1 + n * N)) / 2


def horner(K):
    # symbols = list(sorted(K.free_symbols))
    # k = expand(expand(K))
    # P = np.array([[
    #     int(list(ex.subs(
    #         symb, 2).as_coefficients_dict().values())[0]).bit_length() - 1
    #     for symb in symbols
    # ] for ex in k.as_coefficients_dict()])
    return K


def factor(k):
    if isinstance(k, int):
        return k
    frac = k.as_numer_denom()
    frac = list(map(horner, frac))
    return frac[0] / frac[1]


def quad():
    N1 = N(-1, -1) / 2
    N2 = N(1, -1) / 2
    N3 = N(1, 1) / 2
    N4 = N(-1, 1) / 2
    return [N1, N2, N3, N4]


def serindipity():
    N1 = N(-1, -1) / 2 - (N(-e, -1) + N(-1, -n)) / 2
    N2 = N(1, -1) / 2 - (N(-e, -1) + N(1, -n)) / 2
    N3 = N(1, 1) / 2 - (N(-e, 1) + N(1, -n)) / 2
    N4 = N(-1, 1) / 2 - (N(-e, 1) + N(-1, -n)) / 2
    N5 = N(-e, -1)
    N6 = N(1, -n)
    N7 = N(-e, 1)
    N8 = N(-1, -n)
    return [N1, N2, N3, N4, N5, N6, N7, N8]


def bubble():
    N9 = N(-e, -1)
    N5 = N(-e, -1) - N9
    N6 = N(1, -n) - N9
    N7 = N(-e, 1) - N9
    N8 = N(-1, -n) - N9
    N1 = N(-1, -1) / 2 - (N(-e, -1) + N(-1, -n)) / 2
    N2 = N(1, -1) / 2 - (N(-e, -1) + N(1, -n)) / 2
    N3 = N(1, 1) / 2 - (N(-e, 1) + N(1, -n)) / 2
    N4 = N(-1, 1) / 2 - (N(-e, 1) + N(-1, -n)) / 2
    return [N1, N2, N3, N4, N5, N6, N7, N8, N9]


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

x0 = [A, B, C, D, E, F, G, H]
y0 = [S, T, U, V, W, X, Y, Z]

# We only do the plane strain case
D = [[l + 2 * u, l, 0], [l, l + 2 * u, 0], [0, 0, u]]


def B(N, x, y, J):
    Nx = N.diff(e) * y.diff(n) - N.diff(n) * y.diff(e)
    Ny = -N.diff(e) * x.diff(n) + N.diff(n) * x.diff(e)
    return np.array([[Nx, 0, Ny], [0, Ny, Nx]]).T / J


def format_code(code,
                convert=True,
                padding=13,
                cols=80,
                initial_offset=4,
                sep=" ...\n    "):
    inner_width = cols - padding
    initial_width = inner_width - initial_offset

    if convert:
        code = sympy2.octave_code(code)[2:-2]

    code = code.replace("{", "[").replace("}", "]").replace(".*", "*").replace(
        "], [", "; ").replace("**", "^").replace(".^", "^")
    matches = re.findall("([dx]?[0-9]+[xy]?)", code)
    code = re.sub("sqrt\(", "&", code)
    code = re.sub("\s", "", code)
    code = re.sub("([dx]?[0-9]+)[xy]?", "?", code)

    code = code[:initial_width] + sep + sep.join([
        code[i:i + inner_width]
        for i in range(initial_width, len(code), inner_width)
    ])
    return code.replace('?', "{}").format(*matches).replace("&", "sqrt(")


def symbolic_equations(nodes, dim=2):
    Ns = shapes[nodes]()
    x = sum([n * x_ for n, x_ in zip(Ns, x0[:nodes])])
    y = sum([n * y_ for n, y_ in zip(Ns, y0[:nodes])])

    J = (sympy.diff(x, e) * sympy.diff(y, n) -
         sympy.diff(x, n) * sympy.diff(y, e)).expand()

    # Init K and F to 0
    K = np.zeros((nodes * dim, nodes * dim), dtype="object")
    F = np.zeros((nodes * dim, 1), dtype="object")
    # Init value for our stress calculations
    stress_components = {1: 1, 2: 3, 3: 6}[dim]
    M = np.zeros((nodes * stress_components, nodes * stress_components),
                 dtype="object")
    P = np.zeros((nodes * stress_components, 1), dtype="object")

    # a list of Our integration points and weightings.
    weights = weight_lookup[nodes]
    X = integration_points[weights]
    W = integration_weights[weights]

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
                    m: Xn
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
        p = Na * J
        P[stress_components * i:stress_components * i +
          stress_components] = sum([
              Wx * Wy * p.subs({
                  e: Xe,
                  n: Xn
              }) for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
          ])

    # We put it in horner form since it is supposed to be one of the best forms
    # for polynomial evaluation.
    K = list(map(factor, np.reshape(K, (1, -1))[0]))
    M = list(map(factor, np.reshape(M, (1, -1))[0]))
    F = np.array(
        list(map(lambda f: list(map(sympy2.polys.polyfuncs.horner, f)),
                 F))).T[0]
    P = np.array(
        list(map(lambda p: list(map(sympy2.polys.polyfuncs.horner, p)),
                 P))).T[0]
    precompute, reduced = cse(np.concatenate((K, M, F, P)))

    # Reshape our square matrices
    P = reduced[len(K) + len(M) + len(F):]
    F = reduced[len(K) + len(M):-len(P)]
    M = np.reshape(reduced[len(K):-len(F) - len(P)],
                   (nodes * stress_components, nodes * stress_components))
    K = np.reshape(reduced[:len(K)], (nodes * dim, nodes * dim))

    # Generate the CSE components.
    code = ";".join([
        sympy2.octave_code(v) + "=" + sympy2.octave_code(k)
        for (v, k) in precompute
    ])
    return """
    {}
    K=[{}];
    F=[{}];
    M=[{}];
    P=[{}];""".format(format_code(code, convert=False), format_code(K),
                      format_code(F), format_code(M), format_code(P))


def main():
    renderer = pystache.Renderer()
    quad = symbolic_equations(4)
    serindipity = symbolic_equations(8)
    # Generate our stiffness equations.
    with open('matlab/local_values.m', 'w') as file:
        print(renderer.render_path("templates/local_values.m.tmpl", {
            "quad": quad,
            "serindipity": serindipity
        }),
              file=file)


if __name__ == "__main__":
    main()
