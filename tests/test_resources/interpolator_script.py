import dolfin as df
import pantarei as pr
import numpy as np
import matplotlib.pyplot as plt


from gmri2fem.models.multidiffusion_model import total_concentration

domain = pr.MMSDomain(200)
el = df.FiniteElement("CG", domain.ufl_cell(), degree=2)
element = df.VectorElement(el.family(), el.cell(), el.degree(), dim=2)
W = df.FunctionSpace(domain, element)
V = df.FunctionSpace(domain, el)

exp1 = df.Expression("+(x[0]*x[0] + x[1]*x[1])", degree=2)
exp2 = df.Expression("-(x[0]*x[0] + x[1]*x[1])", degree=2)

compartments = [1, 2]
exps = {1: exp1, 2: exp2}
phi = {1: 0.5, 2: 0.5}
u = pr.assign_mixed_function(exps, W, compartments)
total = total_concentration(W, phi, compartments)
uT = total(u)
assert np.allclose(uT.vector().vec().array, 0.0)