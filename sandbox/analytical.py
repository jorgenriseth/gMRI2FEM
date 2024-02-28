

from scipy.special import erfc 
from scipy import sqrt, exp  

# units: SI:  m, s, kg 
# 3k Dex 
D = 1.03e-10
# 100 um 
z = 100 * 10**-6
# 30 min 
t = 60 * 30
# tortuosity 
lmbda = 1.7 
Deff = D/(lmbda*lmbda)

print ("use units of m, s")  

def robin_analytical(De, z, t, beta, c0 = 1): 
    zeta = beta  
    c = c0 * (erfc(z/(2*sqrt(De*t))) - exp (zeta*z + zeta**2 * De * t) * erfc( z/(2*sqrt(De*t)) + zeta*sqrt(De*t)  ))
    return c 
 

alphas = [0.01, 0.03, 0.1, 0.3,  1, 3, 10, 100] 
for alpha in alphas: 
    dirichlet = erfc(z/(2*sqrt(alpha*Deff*t)))
    print ("when using Dirichlet and ", alpha, " De ", dirichlet, " at ", z)  

print ("Hence at 100 um alpha De ensures significant amounts of tracer unless alpha << 1 ") 

print (""*2)
z = 1000 * 10**-6
alphas = [0.01, 0.03, 0.1, 0.3,  1, 3, 10, 100] 
for alpha in alphas: 
    dirichlet = erfc(z/(2*sqrt(alpha*Deff*t)))
    print ("when using Dirichlet and ", alpha, " De ", dirichlet, " at ", z)  
print ("At 1000 um results change dramatically around alpha=1") 



# z back to start 
z = 100 * 10**-6
print (""*2)
betas = [10**i for i in range(-5, 6)]
for beta in betas:  
    robin = robin_analytical (Deff, z, t, beta) 
    print ("when using Robin with beta ", beta, " ", robin)  
print ("At 100 um, beta less than 10000 really has effect")  

print (""*2)

#betas = [300*i for i in range(1, 100)]
betas = [300, 3000, 30000]
for beta in betas:  
    robin = robin_analytical (Deff, z, t, beta) 
    print ("when using Robin with beta ", beta, " ", robin)  
print ("At 1000 um, beta in (300, 30000) is interesting")  



