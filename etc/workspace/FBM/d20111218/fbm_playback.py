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

# correct code candidate
# same value about date's impliment
def _energy3(v,L):
    term = 0.0
    rowsz, colsz = len(L), len(L[0])
    for i in xrange(rowsz):
        for j in xrange(i,colsz):
            term += v[i]*L[i,j]*v[j]
    return -term

# match one's own way of thinking
def _energy9(v,L):
    term = 0.0
    rowsz, colsz = len(L), len(L[0])
    for i in xrange(rowsz):
        term += L[i,i]
        for j in xrange(i+1,colsz):
            term += v[i]*L[i,j]*v[j]

    return -term

def energy(v,L):
    return _energy9(v,L)
    #return _energy3(v,L)

    #term = dot(dot(v,L),v);
    #return -0.5*term;

def prop(i,v,L,temperature):
    Li = L[i,:]
    term = dot(v,Li);
    return sigma(term/temperature);

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
prefix = sys.argv[1] if sys.argv[1] else "L3x3"
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
