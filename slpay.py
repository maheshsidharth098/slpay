from scipy.optimize import linprog
import numpy as np

s = [10000.0, 20000.0, 30000.0]
r = [0.05, 0.04, 0.03]
N = 5
K = 12
pmax = 1400.0  
# Don't set pmax and N to very small values, otherwise there could be no solution
#=================================================================================

ns = len(s)
c   = np.zeros(ns*N) + K
Aub = np.zeros([ns+N, ns*N])
bub = [-i for i in s] + [pmax for i in range(N)]
xb  = [ (0, pmax) for i in range(ns*N) ]

for i in range(ns):
    for j in range(N):
        ix = i*N + j
        for k in range(1,K+1):
            Aub[i, ix] -= 1.0/((1+r[i]/K)**(j*K+k))
        Aub[j+ns, ix]  = 1.0

res = linprog(c, A_ub=Aub, b_ub=bub,  
              bounds=xb, options={"disp": True})

print res.x.reshape([ns, N])

