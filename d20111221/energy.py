#!/usr/bin/env python
from pylab import *

def sigma(x):
    return 1./(1.+exp(-x))

def energy0(v,L):
    """
    dummy function
    """
    if len(v) == 3:
        return -0.5*(dot(dot(v,L),v) + dot(diag(L),v))
#        val3 = 0.0
#        return -0.5*(dot(dot(v,L),v) - 2.*dot(diag(L),v))
#        return -0.5*(dot(dot(v,L),v) - dot(diag(L),v))
#        return -0.5*(dot(dot(v,L),v)) - trace(L)
#        return -0.5*(dot(dot(v,L),v) + 2.*trace(L))
#         for i in range(0,len(L)):
#             for j in range(i,len(L[0])):
#                 val3 += v[i]*L[i,j]*v[j]
#                 pass
#         return -val3
    
    elif len(v) == 4:
        return -(L[0,1]*v[1] +
                 L[0,2]*v[2] +
                 L[0,3]*v[3] +
                 L[1,2]*v[1]*v[2] +
                 L[1,3]*v[1]*v[3] +
                 L[2,3]*v[2]*v[3] )


def energy1(v,L):
#  using namespace boost::numeric::ublas;
#  double term = inner_prod(prod(v,L),v);
    return -0.5*dot(dot(v,L),v)

def energy2(v,L):
    """
    if \all v is on bit
    0.5*sum(L)
    """
    val2 = 0.0
    for i in range(0,len(L)):
        for j in range(0,len(L[0])):
            val2 += v[i]*L[i,j]*v[j]
            pass
        pass
    return -0.5*val2

def energy3(v,L):
    """
    if \all v is on bit
    upper(L) + diag(L)
    """
    val3 = 0.0
    print L
    assert False
    for i in range(0,len(L)):
        for j in range(i,len(L[0])):
            val3 += v[i]*L[i,j]*v[j]
            pass
    return -val3

def energy4(v,L):
    """
    if \all v is on bit
    upper(L)
    """
    val4 = 0.0
    for i in range(0,len(L)):
        for j in range(i+1,len(L[0])):
            val4 += v[i]*L[i,j]*v[j]
    return -val4


def energy5(v,L):
    """
    if \all v is on bit
    0.5*(sum(L) - diag(L))
    """
    val5 = 0.0
    for i in range(0,len(L)):
        for j in range(0,len(L[0])):
            val5 += v[i]*L[i,j]*v[j]
            pass
        val5 -= v[i]*L[i,i]*v[i]
    pass
    return -0.5*val5

def energy6(v,L):
    """
    if \all v is on bit
    0.5*(sum(L) - diag(L))
    """
    val6 = 0.0
    val6 = dot(dot(v,L),v)
    val6 -= trace(L)
    return -0.5*val6


def energy7(v,L):
    """
    loop version energy6
    """
    val7 = 0.0
    for i in range(0,len(L)):
        for j in range(0,len(L[0])):
            val7 += v[i]*L[i,j]*v[j]
            pass
        val7 -= L[i,i]
        pass
  #val7 -= v[i]*L(i,j)*v[j]
    return -0.5*val7

def denergy0(i,s,W):
    """
    date's resume
    term = \sum_{j=1}^{n} W_{ij} s_{j}
    bias = W_{0i}
    """
    term = 0.0
    bias = W[0,i]
    for j in range(1,len(s)):
        term += W[i,j]*s[j]
    return term + bias

def denergy1(i,s,W):
    """
    simple inner prod
    term = \sum_{j} W_{ij} s_{j}
    """
    return dot(s,W[i,:])

def denergy2(i,s,W):
    """
    simple inner prod (loop version)
    term = \sum_{j
    """
    term = 0.0
    for j in range(len(s)):
        term += W[i,j]*s[j]
    return term

def denergy3(i,s,W):
    """
    duplicate i'th element
    """
#    term, bias = 0.0, W[k,k]
    term, bias = 0.0, W[i,i]
    for j in range(len(s)):
        term += W[i,j]*s[j]
    return term + bias

def denergy4(i,s,W):
    """
    term = \sum_{j\i} W_{i,j} s_{j}
    bias = W_{i,i}
    """
    term, bias = 0.0, 0.0
    for j in range(i):
        term += W[i,j]*s[j]
    for j in range(i+1,len(s)):
        term += W[i,j]*s[j]
    bias = W[i,i]
    return term + bias


def dtob(n,size):
    ret = zeros(size)
    for i in xrange(size):
        ret[size-1-i] = (n>>i) & 0x01
    return ret

def btod(x):
    ret = 0.0
    for i in xrange(size):
        ret += (0x01<<i)*x[size-1-i]
    return ret

D = 3
K=1000
N=1000
alpha = 1.0;

prob_Q = array([0.10, 0.10, 0.05, 0.05, 0.10, 0.10, 0.40, 0.10])
accu_Q = add.accumulate(prob_Q)
all_kurata = [dtob(i,D) for i in xrange(0x01<<D)]
all_date = [dtob(i,D+1) for i in xrange(0x01<<D)]
for ptn in all_date:
    ptn[0] = 1


print "kurata's data list"
for ptn in all_kurata:
    print ptn

print "date's data list"
for ptn in all_date:
    print ptn


L_kurata = zeros((D,D),dtype=float)
L_date = zeros((D+1,D+1),dtype=float)

for n in xrange(N):
    for i,aq in enumerate(accu_Q):
        if n < N*aq:
            break

    L_kurata += outer(all_kurata[i],all_kurata[i])
    L_date += outer(all_date[i],all_date[i])

L_date = L_date - diag(diag(L_date))

print "L_kurata="
print L_kurata

print "L_date="
print L_date

energy_funcs  = [energy0,
                 energy1,
                 energy2,
                 energy3,
                 energy4,
                 energy5,
                 energy6,
                 energy7]

denergy_funcs = [denergy0,
                 denergy1,
                 denergy2,
                 denergy3,
                 denergy4]

for i,energy in enumerate(energy_funcs):
    print "energy%d"%(i)
    print "kurata\t\tdate"
    for v_kurata, v_date in zip(all_kurata,all_date):
        e_kurata = energy(v_kurata, L_kurata)
        e_date = energy(v_date, L_date)
        print "%f\t%f"%(e_kurata,e_date)



for i,denergy in enumerate(denergy_funcs):
    print "denergy%d"%(i)
    print "kurata\t\tdate"
    all_pattern = zip(all_kurata,all_date)
    for n, v in enumerate(all_pattern):
        v_kurata, v_date = v[0],v[1]
        print "n=%d\tkurata"%(n),v_kurata
        print "n=%d\tdate  "%(n),v_date
        for k in xrange(D):
            e_kurata = denergy(k  , v_kurata, L_kurata)
            e_date   = denergy(k+1, v_date  , L_date)
            print "k=%d\t%f\t%f"%(k,e_kurata,e_date)

        print "-------------------------------------"
