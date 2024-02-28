
dt = 0.1 
T = 60*2 

t = 0 
u_prev = 0  
v_prev = 1 
beta = 1e-2 
similar = False 
while not similar: 
    t += dt 
    u = dt * beta * ( v_prev - u_prev) + u_prev 
    v = dt * beta * ( u_prev - v_prev) + v_prev 
    u_prev = u 
    v_prev = v 

    print (t, u, v, u+v) 
    if u > 0.9*v: similar = True 

print ("Final time ", t) 
print ("#"*30) 

rho_u = 0.14 
rho_v = 0.01 
u_prev = 0  
v_prev = 1 
t = 0 
similar = False 
while not similar: 
    t += dt 
    u = (dt / rho_u)  * beta * ( v_prev - u_prev) + u_prev 
    v = (dt / rho_v) * beta * ( u_prev - v_prev) + v_prev 
    u_prev = u 
    v_prev = v 

    print (t, u, v, rho_u*u+rho_v*v) 
    if u > 0.9*v: similar = True 


print ("Final time ", t) 
print ("#"*30) 
rho_u = 0.14 
rho_v = 0.01 
u_prev = 1  
v_prev = 0 
t = 0 
similar = False 
while not similar: 
    t += dt 
    u = (dt / rho_u)  * beta * ( v_prev - u_prev) + u_prev 
    v = (dt / rho_v) * beta * ( u_prev - v_prev) + v_prev 
    u_prev = u 
    v_prev = v 

    print (t, u, v, rho_u*u+rho_v*v) 
    if v > 0.9*u: similar = True 



print ("Final time ", t) 
