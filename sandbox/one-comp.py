
from dolfin import *
import string
import scipy
from scipy.special import erfc 
from numpy import sqrt 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--L_max',default=0.003, type=float)      
parser.add_argument('--Ns',default=[20, 40], type=int, nargs='+')      
parser.add_argument('--dt',default=10, type=int)      
parser.add_argument('--final_time', default=30*60, type=float)
parser.add_argument('--D', default=1.03e-10, type=float)
parser.add_argument('--beta', default=6000, type=float)
parser.add_argument('--alpha', default=1, type=float)
parser.add_argument('--bc', default="Robin", type=str)
parser.add_argument('--check_points', default=[0.0001, 0.001], type=float, nargs='+')
parser.add_argument('--lmbda', default=1.7, type=float)
args = parser.parse_args()

print ("args ", args)

numerical_solutions = []
for N in args.Ns: 
    print ("solving for ", N)

    mesh = UnitIntervalMesh(N)
    mesh.coordinates()[:] *= args.L_max  
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V) 


    # units are m and seconds
    Deff = Constant(args.D/(args.lmbda*args.lmbda))
    dt = Constant(args.dt)
    beta = Constant(args.beta) 
    alpha = Constant(args.alpha) 

    
    # assume these are set to zero
    U = Function(V)
    U_prev = Function(V)
    bc_func = Constant(1)  

    class Left(SubDomain): 
      def inside(self, x, on_boundary): 
        return near(x[0], 0.0)
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    left = Left()
    left.mark(boundaries, 2)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    a =  u*v*dx + dt*Deff*inner(grad(u), grad(v))*dx 
    L = dt*Constant(0)*v*dx + U_prev*v*dx  
    if args.bc == "Robin": 
      print ("Robin bc") 
      L += dt*args.beta*Deff*(bc_func*v*ds(2)) 
      a += dt*args.beta*Deff*(u*v*ds(2))

    MA = assemble(a) 

    t = 0.0
    bc = 0 
    if args.bc == "Dirichlet": 
      bc = DirichletBC(V, bc_func, left_boundary)
      bc.apply(MA)


    ts = []
    u0s = []
    uls = []
    urs = []
    done = False
    while t < args.final_time: 
        bc_func.t = t 
        t += args.dt  
        ts.append(t)
        b = assemble(L) 
        if args.bc == "Dirichlet": bc.apply(b) 
        solve(MA,U.vector(), b, "gmres", "ilu")
        U_prev.assign(U)
    numerical_solutions.append((V.tabulate_dof_coordinates(), U.vector()))

    for check_point in args.check_points: 
      v = U(check_point) 
      print ("tracer at check_point ", check_point, v, " at time ", t, " alpha ", args.alpha, " beta ", args.beta)


x = V.tabulate_dof_coordinates()
xx = scipy.arange(0, args.L_max, args.L_max / 1000) 


Deffv = Deff.values()[0]
analytical_solution = erfc(xx / (2*sqrt( Deffv * t))) 
if args.bc == "Robin": 
    from scipy.special import erfc 
    from scipy import sqrt, exp  

    zeta = args.beta  
    c0 = 1.0 
    z = xx
    c = c0 * (erfc(z/(2*sqrt(Deffv*t))) - exp (zeta*z + zeta**2 * Deffv * t) * erfc( z/(2*sqrt(Deffv*t)) + zeta*sqrt(Deffv*t)  ))
    analytical_solution = c 


import matplotlib.pyplot as plt 
plt.gcf().subplots_adjust(bottom=0.25, left=0.25)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.xlabel("distance [mm]", size=24)
plt.ylabel("concentration", size=24)

for x, vec in numerical_solutions: 
    plt.plot(x, vec, linewidth=7)
plt.plot(xx, analytical_solution, linewidth=7)
legend = ["N=%d"%N for N in args.Ns]
legend.append("analytical sol")
plt.legend(legend,  prop={"size" : 20}, loc=1)
file_str = "3kDex_numerical_1D_L_max%d_dt%d_final%d"% (args.L_max,args.dt, args.final_time)
plt.savefig(file_str)
plt.show()





