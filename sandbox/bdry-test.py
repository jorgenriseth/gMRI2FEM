# import time as pytime

import dolfin as df
import matplotlib.pyplot as plt
import numpy as np
import pantarei as pr
from ufl import dot, ds, dx, grad, inner, lhs, rhs

domain = df.UnitIntervalMesh(100)
V = df.FunctionSpace(domain, "CG", 1)

dt = 0.001
time = pr.TimeKeeper(dt=dt, endtime=1)

U = df.Expression(f"{np.pi}*sin(2*{np.pi}*t)", domain=domain, degree=2, t=time)
c_expr = df.Expression(
    f"exp(-pow(((x[0] - (1 - cos(2*{np.pi}*t))/ 2))/l, 2))",
    l=0.5,
    # domain=domain,
    degree=2,
    t=time,
)
n = df.Expression(f"x[0] <= 0.5 ? -1 : 1", degree=0)


def left_boundary(x, on_boundary):
    if on_boundary:
        if df.near(x[0], 0, 1e-14):
            return True
    return False


def right_boundary(x, on_boundary):
    if on_boundary:
        if df.near(x[0], 1, 1e-14):
            return True
    return False


c, w = df.TrialFunction(V), df.TestFunction(V)
c0 = df.interpolate(c_expr, V)
F = ((c - c0) / dt) * w * dx - inner(c * U, w.dx(0)) * dx  # - w * c * U * n * ds
a = lhs(F)
l = rhs(F)

bcs = [
    # df.DirichletBC(V, c_expr, left_boundary),
    df.DirichletBC(V, c_expr, "on_boundary"),
]

c = df.Function(V)
for idx, ti in enumerate(time):
    # plt.figure()
    print(float(time))
    A = df.assemble(a)
    b = df.assemble(l)
    for bc in bcs:
        bc.apply(A, b)
    df.solve(A, c.vector(), b)
    if idx % 10 == 0:
        df.plot(c)
        df.plot(c_expr, mesh=domain)
    plt.ylim(0, 1.2)
    plt.title(f"t={float(time)}")
    c0.assign(c)
plt.show()
