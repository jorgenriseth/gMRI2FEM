
from dolfin import *
import string
import scipy
from scipy.special import erfc 
from numpy import sqrt 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--L_max',default=0.003, type=float)      
parser.add_argument('--Ns',default=[80], type=int, nargs='+')      
parser.add_argument('--dt',default=10, type=float)      
parser.add_argument('--final_time', default=30*60, type=float)
parser.add_argument('--D', default=1.03e-10, type=float)
parser.add_argument('--gamma', default=1e-2, type=float)
parser.add_argument('--beta1', default=3000, type=float)
parser.add_argument('--beta2', default=300000, type=float)
parser.add_argument('--alpha', default=2, type=float)
parser.add_argument('--bc', default="Robin", type=str)
parser.add_argument('--check_points', default=[0.0001, 0.001], type=float, nargs='+')
parser.add_argument('--lmbda', default=1.7, type=float)
parser.add_argument('--rho1', default=0.14, type=float)
parser.add_argument('--rho2', default=0.01, type=float)
parser.add_argument('--v', default=2.0e-5, type=float)
parser.add_argument('--sigma', default=5.0e-5, type=float) 

args = parser.parse_args()

print ("args ", args)

numerical_solutions2 = []
for N in args.Ns: 
    print ("solving for ", N)
    print ("U1 is ECS while U2 is PVS") 

    mesh = UnitIntervalMesh(N)
    mesh.coordinates()[:] *= args.L_max  
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    P12 = P1*P1  
    W = FunctionSpace(mesh, P12)
    V = FunctionSpace(mesh, P1) 
    (u1,u2) = TrialFunctions(W)
    (v1,v2) = TestFunctions(W) 


    D = Constant(args.D)
    D1eff = Constant(args.D/(args.lmbda*args.lmbda))
    D2eff = Constant(args.D*args.alpha)
    dt = Constant(args.dt)
    beta1 = Constant(args.beta1) 
    beta2 = Constant(args.beta2) 
    gamma = Constant(args.gamma) 
    rho1 = Constant(args.rho1) 
    rho2 = Constant(args.rho2) 
    sigma = Constant(args.sigma)
    v = Constant(args.v) 

    class Left(SubDomain): 
      def inside(self, x, on_boundary): 
        return near(x[0], 0.0)


    bc_func1 = Constant(0)  
    bc_func2 = Constant(1)  
#    bcs = [DirichletBC(W.sub(0), bc_func, left_boundary), DirichletBC(W.sub(1), bc_func, left_boundary)]
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    left = Left()
    left.mark(boundaries, 2)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)


    # assume these are set to zero
    U = Function(W)
    U1 = Function(V) 
    U2 = Function(V) 
    U_tot = Function(V)
    U_prev= Function(W)
    U1_prev = Function(V) 
    U2_prev = Function(V) 

    a11 =  u1*v1*dx + dt*D1eff*inner(grad(u1), grad(v1))*dx 
    a11 += dt*beta1*D1eff*(u1*v1*ds(2))
    a11 += dt*gamma*(u1-u2)*v1*dx 
    a22 =  u2*v2*dx + dt*D2eff*inner(grad(u2), grad(v2))*dx 
    a22 += dt*beta2*D2eff*(u2*v2*ds(2))
    a22 += dt*gamma*(u2-u1)*v2*dx 
    a22 += dt*sigma*u2*v2*dx 
    a22 += dt*v*u2.dx(0)*v2*dx  
    a  =  a11 + a22  
 
    L1 = dt*Constant(0)*v1*dx + U1_prev*v1*dx  
    L1 += dt*beta1*D1eff*(bc_func1*v1*ds(2)) 

    L2 = dt*Constant(0)*v2*dx + U2_prev*v2*dx  
    L1 += dt*beta2*D2eff*(bc_func2*v2*ds(2)) 
    L = L1 + L2 

    MA = assemble(a) 


    ts = []
    u0s = []
    uls = []
    urs = []
    done = False
    t = 0.0
    while t < args.final_time: 
        bc_func1.t = t 
        bc_func2.t = t 
        t += args.dt  
        ts.append(t)
        b = assemble(L) 

        solve(MA, U.vector(), b, "gmres", "ilu")
        U_prev.assign(U)
        assigner = FunctionAssigner([V, V], W)
        assigner.assign([U1_prev, U2_prev], U_prev)  
        assigner.assign([U1, U2], U)  
        U_tot.vector()[:] = rho1.values()[0]*U1.vector()[:] + rho2.values()[0]*U2.vector()[:] 


    numerical_solutions2.append((V.tabulate_dof_coordinates(), U1.vector()))
    numerical_solutions2.append((V.tabulate_dof_coordinates(), U2.vector()))
    numerical_solutions2.append((V.tabulate_dof_coordinates(), U_tot.vector()))

    for check_point in args.check_points: 

      v = U1(check_point) 
      print ("U1   tracer at check_point ", check_point, v,   " at time ", t)
      v = U2(check_point) 
      print ("U2   tracer at check_point ", check_point, v,   " at time ", t)
      v = U_tot(check_point) 
      print ("Utot tracer at check_point ", check_point, v,   " at time ", t)






import matplotlib.pyplot as plt 
plt.gcf().subplots_adjust(bottom=0.25, left=0.25)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.xlabel("distance [mm]", size=24)
plt.ylabel("concentration", size=24)


for x, vec in numerical_solutions2: 
    plt.plot(x, vec, linewidth=7)
legend = [("U1 N=%d"%N, "U2 N=%d"%N, "Ut N=%d"%N)   for N in args.Ns]
legend = [item for sublist in legend for item in sublist]
plt.legend(legend,  prop={"size" : 20}, loc=1)
file_str = "Two_comp_3kDex_numerical_1D_L_max%d_dt%d_final%d"% (args.L_max,args.dt, args.final_time)
plt.savefig(file_str)
plt.show()





