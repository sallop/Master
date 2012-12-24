from pylab import *

def dtob(n, size):
    ret = zeros(size,dtype=float)
    for i in xrange(size):
        ret[size-1-i] = (n>>i) & 0x01;

    return ret;


def btod(x):
    ret  = 0
    size = len(x)
    for i in xrange(size):
        ret += x[size-1-i]*(0x01 << i);
    return ret;


D = 3
N = 1000
gamma = [dtob(i,D) for i in xrange(0x01<<D)]

prob_Q = [0.10, 0.10, 0.05, 0.05, 0.10, 0.10, 0.40, 0.10]
accm_Q = add.accumulate(prob_Q)
accm_Q[-1]=1.0

L=zeros((D,D),dtype=float)

for n in xrange(N):
    #r = rand()
    #for k,aq in enumerate(accm_Q):
    #    if r < aq: break
    for k,aq in enumerate(accm_Q):
        if n < N*aq: break

    L += outer(gamma[k],gamma[k])
    

L /= N

# result
# array([[ 0.7 ,  0.5 ,  0.2 ],
#        [ 0.5 ,  0.6 ,  0.15],
#        [ 0.2 ,  0.15,  0.35]])

# In [3]: linalg.eig(L)
# Out[3]: 
# (array([ 1.22348167,  0.1438558 ,  0.28266252]),
#  array([[ 0.716947  ,  0.69088651, -0.09307427],
#        [ 0.6409281 , -0.70576658, -0.30183556],
#        [ 0.27422283, -0.15674618,  0.94880581]]))
