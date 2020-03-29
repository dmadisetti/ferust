import sympy
from sympy.abc import n, e, h, l, u, h, A, B, C, D, E, F, G, H, S, T, U, V, W, X, Y, Z
from sympy.printing import octave_code
import numpy as np
import pystache
from itertools import product
import re

# using the shape function convention in #5
N = lambda E, N: ((1 + e * E) * (1 + n * N)) / 2


def quad():
    N1 = N(-1, -1) / 2
    N2 = N(1, -1) / 2
    N3 = N(1, 1) / 2
    N4 = N(-1, 1) / 2
    return [N1, N2, N3, N4]


def serindipity():
    N1 = sympy.factor(
        sympy.simplify(N(-1, -1) / 2 - (N(-e, -1) + N(-1, -n)) / 2))
    N2 = sympy.factor(
        sympy.simplify(N(1, -1) / 2 - (N(-e, -1) + N(1, -n)) / 2))
    N3 = sympy.factor(sympy.simplify(N(1, 1) / 2 - (N(-e, 1) + N(1, -n)) / 2))
    N4 = sympy.factor(
        sympy.simplify(N(-1, 1) / 2 - (N(-e, 1) + N(-1, -n)) / 2))
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
    combination(-sympy.sqrt(3) / sympy.sqrt(5), 0,
                sympy.sqrt(3) / sympy.sqrt(5))
}
integration_weights = {
    4:
    combination(1, 1),
    8:
    combination(
        sympy.Rational(5) / 9,
        sympy.Rational(8) / 9,
        sympy.Rational(5) / 9)
}

x0 = [A, B, C, D, E, F, G, H]
y0 = [S, T, U, V, W, X, Y, Z]

# We only do the plane strain case
lb = l  #2 * l * u / (l + 2 * u)
D = [[lb + 2 * u, lb, 0], [lb, lb + 2 * u, 0], [0, 0, u]]


def B(N):
    return np.array([[N.diff(e), 0, N.diff(n)], [0, N.diff(n), N.diff(e)]]).T


def format_code(K, padding=13, cols=80, initial_offset=4, sep=" ...\n    "):
    inner_width = cols - padding
    initial_width = inner_width - initial_offset

    code = octave_code(K)[2:-2].replace("{", "[").replace("}", "]").replace(
        ".*", "*").replace("], [", "; ")
    c = code
    matches = re.findall("[^0-9|^]([0-9]+)[^0-9|$]", code)
    code = re.sub("([0-9]+)", "?", code)

    code = code[:initial_width] + sep + sep.join([
        code[i:i + inner_width]
        for i in range(initial_width, len(code), inner_width)
    ])
    print(code)
    print("\n".join(matches))
    return code.replace('?', "{}").format(*matches)


def symbolic_stiffness(nodes, dim=2):
    Ns = shapes[nodes]()
    x = sum([n * x_ for n, x_ in zip(Ns, x0[:nodes])])
    y = sum([n * y_ for n, y_ in zip(Ns, y0[:nodes])])

    J = sympy.diff(x, e) * sympy.diff(y, n) - sympy.diff(x, n) * sympy.diff(
        y, e)

    # Init K to 0
    K = np.zeros((nodes * dim, nodes * dim), dtype="object")

    # a list of Our integration points and weightings.
    X = integration_points[nodes]
    W = integration_weights[nodes]

    for i, Na in enumerate(Ns):
        for j, Nb in enumerate(Ns):
            # symmetric
            if i > j:
                continue
            k = sympy.simplify(
                np.reshape(
                    sympy.simplify(np.matmul(np.matmul(B(Na).T, D), B(Nb))),
                    (dim, dim)))
            k = np.reshape(
                np.sum([
                    Wx * Wy * k.subs({
                        e: Xe,
                        n: Xn
                    }) * J.subs({
                        e: Xe,
                        n: Xn
                    }) for ((Xe, Xn), (Wx, Wy)) in zip(X, W)
                ],
                       axis=0), (dim, dim))
            K[2 * i:2 * i + 2, 2 * j:2 * j + 2] += k

            # So we end up with just URH
            # NOTE: Only works in 2D
            if i == j:
                K[i + 1, i] = 0

    # K = list(map(lambda k: list(map(sympy.polys.polyfuncs.horner, k)), K))
    K = sympy.simplify(K)
    return format_code(K)


# code = symbolic_stiffness(4)

renderer = pystache.Renderer()
print(
    renderer.render_path("local_stiffness.m.tmpl", {
        "quad": symbolic_stiffness(4),
        "serindipity": symbolic_stiffness(8)
    }))
