#!/usr/bin/env python
#-*- coding: utf-8 -*-
from pylab import *
import viewer as bv
import matplotlib.pyplot as plt
import hinton_dia1

from fbm_option import *
from readcfg import *
from read_ublas import *

def sigma(x):
    return 1./(1+exp(-x))

def _energy0(v,L):
    ene=0.0
    for i in xrange(1,len(L)):
        ene += v[i]*L[0,i]
        for j in xrange(i+1,len(L[0])):
            ene += v[i]*L[i,j]*v[j]
    return -ene

def _energy1(v,L):
    term = dot(dot(v,L),v)
    return -0.5*term

def _energy2(v,L):
    val2 = 0.0
    for i in xrange(0,len(L)):
        for j in xrange(j,len(L[0])):
            val2 += v[i]*L[i,j]*v[j]
    return -0.5*val2

# correct code candidate
def _energy3(v,L):
    val3 = 0.0
    for i in xrange(0,len(L)):
        for j in xrange(i,len(L[0])):
            val3 += v[i]*L[i,j]*v[j]
    return -val3

def _energy4(v,L):
    val4 = 0.0
    for i in xrange(0,len(L)):
        for j in xrange(i+1,len(L[0])):
            val4 += v[i]*L[i,j]*v[j]
    return -val4


def _energy5(v,L):
    val5 = 0.0
    for i in xrange(0,len(L)):
        val5 -= v[i]*L[i,i]*v[i]
        for j in xrange(0,len(L[0])):
            val5 += v[i]*L[i,j]*v[j]
    return -0.5*val5

def _energy6(v,L):
    val6 = dot(dot(v,L),v) - diag(diag(L))
    return -0.5*val6

def _energy7(v,L):
    val7 = 0.0
    for i in xrange(0,len(L)):
        val7 -= L[i,i]
        for j in xrange(0,len(L[0])):
            val7 += v[i]*L[i,j]*v[j]
    return -0.5*val7

# correct code candidate
# same value at date's impliment
def _energy8(v,L):
    ret = 0.0
    for i in xrange(0,len(L)):
        for j in xrange(i,len(L[0])):
            ret += v[i]*L[i,j]*v[j]
    return -ret

# I think this is hold water.
# Kurata's Report 「ボルツマン・マシン」 describe this
def _energy9(v,L):
    ret=0.0
    for i in xrange(0,len(L)):
        for j in xrange(i+1,len(L[0])):
            ret += v[i]*L[i,j]*v[j]
    ret += L[i,i]
    return -ret

# I think this is hold water 2. 
# 「ニューロンコンピューティングの基礎理論」
def _energy10(v,L):
    ret=0.0
    for i in xrange(0,len(L)):
        ret += v[i]*L[i,i]
        for j in xrange(i+1,len(L[0])):
            ret += v[i]*L[i,j]*v[j]
    return -ret

def _denergy0(i,v,L):
    u_i=0.0
    for j in xrange(1,len(v)):
        u_i += L[i,j]*v[j]
    return u_i + L[i,i]

def _denergy1(i,v,L):
    u_i=0.0
    for j in xrange(0,len(v)):
        u_i += L[i,j]*v[j]
    return u_i - L[i,i]*v[i] + L[i,i]

# correct code candidate
# same value date's impliments
def _denergy2(i,v,L):
    u_i=0.0
    for j in xrange(0,len(v)):
        u_i += L[i,j]*v[j]
    return u_i - L[i,i]*v[i]

def _denergy3(i,v,L):
    u_i=0.0
    for j in xrange(0,len(v)):
        u_i += L[i,j]*v[j]
    return u_i

# I think this is hold water,
# at Kurata's model
def _denergy4(i,v,L):
    u_i=0.0
    for j in xrange(0,i):
        u_i += L[i,j]*v[j]
    for j in xrange(i+1,len(v)):
        u_i += L[i,j]*v[j]
    u_i += L[i,i];# bias term
    return u_i

def _prop0(i,v,L,temperature):
    return sigma(_denergy0(i,v,L)/temperature)

def _prop1(i,v,L,temperature):
    return sigma(_denergy1(i,v,L)/temperature)

def _prop2(i,v,L,temperature):
    return sigma(_denergy2(i,v,L)/temperature)

def _prop3(i,v,L,temperature):
    return sigma(_denergy3(i,v,L)/temperature)

def _prop4(i,v,L,temperature):
    return sigma(_denergy4(i,v,L)/temperature)

g_energy_tbl = [
    _energy0, _energy1, _energy2, _energy3, _energy4,
    _energy5, _energy6, _energy7, _energy8, _energy9,
    _energy10]

g_denergy_tbl = [_denergy0, _denergy1, _denergy2,_denergy3,_denergy4,]
g_prop_tbl = [_prop0,_prop1,_prop2,_prop3,_prop4]

g_energy_idx = 0
g_denergy_idx = 0
g_prop_idx = 0

def energy(v,L):
    return g_energy_tbl[g_energy_idx](v,L)

def denergy(i,v,L):
    return g_denergy_tbl[g_denergy_idx](i,v,L)

def prop(i,v,L,temperature):
    return g_prop_tbl[g_prop_idx](i,v,L,temperature)

#// Z(\theta)
def partition_function(L,all_v):
    ret=0.0;
    f = lambda v: exp(-energy(v,L))
    for v in all_v: ret += f(v)
    return ret

def boltzmann_distribution(v,L,all_v):
    numer = exp(-energy(v,L))
    denom = partition_function(L,all_v)
    return numer/denom


def KLd(P,Q):
    d=0.0
    for p,q in zip(P,Q):d += p*log(p/q)
    return d

def print_matrix(W):
    for wr in W:
        print wr

def kstep_gibbs(K,v,L,temperature):
    D = len(v)
    for k in xrange(K):
        i = rand()%D;
        v[i] = 1 if prop(i, v, L, t) > rand() else 0

def free_run_expander():
    return

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

def hinton2(W, maxWeight=None, ax=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    """
    reenable = False
    if plt.isinteractive():
        plt.ioff()

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if not maxWeight:
        maxWeight = 2**ceil(log(abs(W).max())/log(2))

    #plt.clf()
    ax.cla()
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal','box')
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_locator(NullLocator())

    for (x,y),w in ndenumerate(W):
        color = 'white' if w > 0 else 'black'
        size = sqrt(abs(w))
        rect = Rectangle([x-size/2, y-size/2],
                         size,
                         size,
                         facecolor = color,
                         edgecolor = color)
        ax.add_patch(rect)
    ax.autoscale_view()
    ax.set_ylim(*ax.get_ylim()[::-1])
    if reenable:
        plt.ion()

    # Reverse the yaxis limits
    

def hinton_sub(W, maxWeight=None, ax=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if not maxWeight:
        maxWeight = 2**ceil(log(abs(W).max())/log(2))

#    plt.cla()
    ax.cla()
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal','box')
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_locator(NullLocator())

    for (x,y),w in ndenumerate(W):
        color = 'white' if w > 0 else 'black'
        size = sqrt(abs(w))
        rect = Rectangle([x-size/2,y-size/2],
                         size,
                         size,
                         facecolor=color,
                         edgecolor=color)
        ax.add_patch(rect)
    ax.autoscale_view()
    # Reverse the yaxis limits
    #ax.set_ylim(*ax.get_ylim()[::-1])

def hinton_sub2(W, maxWeight=None, ax=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if not maxWeight:
        maxWeight = 2**ceil(log(abs(W).max())/log(2))

    ax.cla()
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal','box')
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_locator(NullLocator())

    for (x,y),w in ndenumerate(W):
        color = 'white' if w > 0 else 'black'
        size = sqrt(abs(w))
        rect = Rectangle([x-size/2,y-size/2],
                         size,
                         size,
                         facecolor=color,
                         edgecolor=color)
        ax.add_patch(rect)
    ax.autoscale_view()
    # Reverse the yaxis limits
    if hinton_sub2.count == 0:
        ax.set_ylim(*ax.get_ylim()[::-1])
        hinton_sub2.count = 1

hinton_sub2.count = 0

dir = "dir-dat"
prefix = sys.argv[1] if sys.argv[1] else "K23x3"
identifies = ["cfg","mat","t_gibbs"]
suffix = ".dat"
files = dict((id, "%s/%s-%s%s"%(dir,prefix,id,suffix))
             for id in identifies)

symcfg = {}
fr = open(files['cfg'],'r')
for i,line in enumerate(fr):
    if i > 19: break
    read_config(line, symcfg)
fr.close()

print "symcfg=",symcfg
for k,v in symcfg.iteritems():
    print k,
    print v

D = int(symcfg['D'].pop())
P = int(symcfg['P'].pop())
K = int(symcfg['K'].pop())
M = int(symcfg['M'].pop())
N = int(symcfg['N'].pop())
T = int(symcfg['T'].pop())
g_energy_idx = int(symcfg['g_energy_idx'].pop())
g_denergy_idx= int(symcfg['g_denergy_idx'].pop())
g_prop_idx   = int(symcfg['g_prop_idx'].pop())

alpha  = float(symcfg['alpha'].pop())
prob_Q = map(lambda v:read_ublas_vec(v), symcfg['prob_Q'  ]).pop()
all_v  = map(lambda v:read_ublas_vec(v), symcfg['all_v[i]'])
prob_P = zeros(0x01<<D)

print "D=",D
print "P=",P
print "K=",K
print "M=",M
print "N=",N
print "T=",T
print "g_energy_idx=",g_energy_idx
print "g_denergy_idx=",g_denergy_idx
print "g_prop_idx=",g_prop_idx
print "alpha=",alpha
print "prob_Q=",prob_Q
print "all_v=",all_v 
#assert False

plt.ion()
fig = plt.figure()
axis = [fig.add_subplot(311),
        fig.add_subplot(312),
        fig.add_subplot(313),]
axis[2].set_xlim(-1,16)
kl_result = []

for line in open(files['mat'],"r"):
    import time

    t, mat = read_t_mat(line)
    L = array(mat, dtype=float)
    prob_P = [boltzmann_distribution(v,L,all_v) for v in all_v]
    ene_result = map(lambda v : energy(v,L), all_v)
    kl_result.append(KLd(prob_P, prob_Q))

    hinton_sub   (L, ax = axis[0])
    bv.cmp_energy(ene_result, ax = axis[1])
    bv.print_KL  (arange(t), kl_result, ax = axis[2])

    plt.draw()

plt.show()
