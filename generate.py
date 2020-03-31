#!/usr/bin/python3
import sympy
from sympy.abc import n, e, h, f, l, u, h, A, B, C, D, E, F, G, H, S, T, U, V, W, X, Y, Z
from sympy.printing import octave_code
import numpy as np
import pystache
from itertools import product
import re
import pdb

# using the shape function convention in #5
N = lambda E, N: ((1 + e * E) * (1 + n * N)) / 2


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


combination = lambda *a: list(product(a, a))
shapes = {4: quad, 8: serindipity}
integration_points = {
    4:
    combination(-1 / sympy.sqrt(3), 1 / sympy.sqrt(3)),
    8:
    combination(-sympy.sqrt(3) / sympy.sqrt(5),
                sympy.sqrt(3) / sympy.sqrt(5), 0)
}
integration_weights = {
    4:
    combination(1, 1),
    8:
    combination(
        sympy.Rational(5) / 9,
        sympy.Rational(5) / 9,
        sympy.Rational(8) / 9)
}

x0 = [A, B, C, D, E, F, G, H]
y0 = [S, T, U, V, W, X, Y, Z]

# We only do the plane strain case
# l = 2 * l * u / (l + 2 * u)
D = [[l + 2 * u, l, 0], [l, l + 2 * u, 0], [0, 0, u]]


def B(N, x, y, J):
    Nx = N.diff(e) * y.diff(n) - N.diff(n) * y.diff(e)
    Ny = -N.diff(e) * x.diff(n) + N.diff(n) * x.diff(e)
    return np.array([[Nx, 0, Ny], [0, Ny, Nx]]).T / J


def format_code(K, padding=13, cols=80, initial_offset=4, sep=" ...\n    "):
    inner_width = cols - padding
    initial_width = inner_width - initial_offset

    code = octave_code(K)[2:-2].replace("{", "[").replace("}", "]").replace(
        ".*", "*").replace("], [", "; ").replace("**", "^")
    matches = re.findall("([0-9]+)", code)
    code = re.sub("sqrt\(", "&", code)
    code = re.sub("([0-9]+)", "?", code)

    code = code[:initial_width] + sep + sep.join([
        code[i:i + inner_width]
        for i in range(initial_width, len(code), inner_width)
    ])
    return code.replace('?', "{}").format(*matches).replace("&", "sqrt(")


def symbolic_equations(nodes, dim=2):
    Ns = shapes[nodes]()
    x = sum([n * x_ for n, x_ in zip(Ns, x0[:nodes])])
    y = sum([n * y_ for n, y_ in zip(Ns, y0[:nodes])])

    J = sympy.diff(x, e) * sympy.diff(y, n) - sympy.diff(x, n) * sympy.diff(
        y, e)

    # Init K and F to 0
    K = np.zeros((nodes * dim, nodes * dim), dtype="object")
    F = np.zeros((nodes * dim, 1), dtype="object")

    # a list of Our integration points and weightings.
    X = integration_points[nodes]
    W = integration_weights[nodes]

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

        # We can concurrently solve for the base F component.
        # We only set the body force from the y component
        temp = Na * f * J
        print(F.shape, 1 + i * 2)
        F[1 + i * 2] = sum([
            Wx * Wy * temp.subs({
                e: Xe,
                n: Xn
            }) for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
        ])

    # We put it in horner form since it is supposed to be one of the best forms
    # for polynomial evaluation.
    K = list(map(lambda k: list(map(sympy.polys.polyfuncs.horner, k)), K))
    F = list(map(lambda f: list(map(sympy.polys.polyfuncs.horner, f)), F))
    return format_code(K), format_code(F)


def main():
    renderer = pystache.Renderer()
    quad_stiffness, quad_force = symbolic_equations(4)
    # Generate our stiffness equations.
    with open('matlab/local_stiffness.m', 'w') as file:
        print(
            renderer.render_path(
                "templates/local_stiffness.m.tmpl",
                {
                    "quad": quad_stiffness,
                    "serindipity": ""  #symbolic_stiffness(8)
                }),
            file=file)
    # Generate our force equations.
    with open('matlab/local_force.m', 'w') as file:
        print(
            renderer.render_path(
                "templates/local_force.m.tmpl",
                {
                    "quad": quad_force,
                    "serindipity": ""  #symbolic_stiffness(8)
                }),
            file=file)


if __name__ == "__main__":
    main()
