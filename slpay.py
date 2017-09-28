from scipy.optimize import linprog
import numpy as np

s  = [10000.0, 40000.0, 5000.0, 5000.0]
r  = [0.055, 0.045, 0.05, 0.06]
nyear  = 5
K      = 1
t0     = 201709
pmax   = 1300.0

# Don't set pmax and nyear to very small values, otherwise there could be no solution
#=================================================================================

N = nyear * 12/K
time = []
for n in range(N):
    if (t0%100+n*K)%12 == 0:
        t = (t0//100 + (t0%100+n*K)//12 - 1)*100 + 12
    else:
        t = (t0//100 + (t0%100+n*K)//12)*100 + (t0%100+n*K)%12
    time.append(t)
    

ns = len(s)
c   = np.zeros(ns*N) + K
Aub = np.zeros([ns+N, ns*N])
bub = [-i for i in s] + [pmax for i in range(N)]
xb  = [ (0, pmax) for i in range(ns*N) ]

for i in range(ns):
    for j in range(N):
        ix = i*N + j
        for k in range(1,K+1):
            Aub[i, ix] -= 1.0/((1+r[i]/12)**(j*K+k))
        Aub[j+ns, ix]  = 1.0

res = linprog(c, A_ub=Aub, b_ub=bub,  
              bounds=xb, options=dict(bland=True, tol=1.0e-6, disp=True))
print 'Maximum monthly payment: ', pmax
print '--------------------------------------'
print ['YYYYMM'] + [' S'+format(i+1) for i in range(ns)]
print '--------------------------------------'
print np.append(np.array(time).reshape([N,1]), 
                res.x.reshape([ns, N]).T, axis=1).astype(int)
