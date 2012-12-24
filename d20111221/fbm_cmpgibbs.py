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

class SimConfig:
    """
    Simulation Configurations
    """
    def __init__(self, prefix):
        self.dir = "dir-dat"
        self.prefix = prefix
        self.identifies = ["cfg","mat","t_gibbs"]
        self.suffix = ".dat"
        self.files = dict((id, "%s/%s-%s%s"%(self.dir, prefix, id, self.suffix))
                          for id in self.identifies)
        self.frs = dict((id, None) for id in self.identifies)
        
        symcfg = {}
        fr = open(self.files['cfg'],'r')
        for i, line in enumerate(fr):
            if i > 32: break
            read_config(line,symcfg)
        fr.close()

        self.D = int(symcfg['D'].pop())
        self.P = int(symcfg['P'].pop())
        self.K = int(symcfg['K'].pop())
        self.M = int(symcfg['M'].pop())
        self.N = int(symcfg['N'].pop())
        self.T = int(symcfg['T'].pop())
        self.g_energy_idx = int(symcfg['g_energy_idx'].pop())
        self.g_denergy_idx= int(symcfg['g_denergy_idx'].pop())
        self.g_prop_idx   = int(symcfg['g_prop_idx'].pop())
        self.alpha  = float(symcfg['alpha'].pop())
        self.prob_Q = map(lambda v:read_ublas_vec(v), symcfg['prob_Q'  ]).pop()
        self.all_v  = map(lambda v:read_ublas_vec(v), symcfg['all_v[i]'])

    def open(self, fname):
        self.frs[fname] = open(self.files[fname],'r')
        return self.frs[fname]

    def readline(self, fname):
        return self.frs[fname].readline()

    def close(self, fname):
        close(self.frs[fname])
        

symconfigs = [SimConfig(prefix) for prefix in sys.argv[1:]]
    
for symcfg in symconfigs:
    print type(symcfg), symcfg
    print "prefix=",symcfg.prefix
    print "D=",symcfg.D
    print "P=",symcfg.P
    print "K=",symcfg.K
    print "M=",symcfg.M
    print "N=",symcfg.N
    print "T=",symcfg.T
    print "g_energy_idx=" ,symcfg.g_energy_idx
    print "g_denergy_idx=",symcfg.g_denergy_idx
    print "g_prop_idx="   ,symcfg.g_prop_idx
    print "alpha="        ,symcfg.alpha
    print "prob_Q="       ,symcfg.prob_Q
    print "all_v="        ,symcfg.all_v 


def printKL(t,kl,ax=None):
    """
    plot KL divergence
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add.subplot(111)

    ax.cla()
    ax.plot(t,kl)
    ax.set_ylim(0.0,1.1)

figsize = (12,8)


#plt.ion()
#plt.ioff()
fig = plt.figure(figsize=figsize)
row, column = len(symconfigs), 2
axis = [fig.add_subplot(row, column, i+1)
        for i in xrange(row*column)]

for symcfg in symconfigs:
    symcfg.open('t_gibbs')

line = ""
kl = dict((symcfg.prefix,[]) for symcfg in symconfigs)


j=0
while True:
#while j < 10:
    import time

    for i,symcfg in enumerate(symconfigs):
        key = symcfg.prefix
        line = symcfg.readline('t_gibbs')
        if not line:
            break
        
        t, gibbs_P = read_t_gibbs(line)
        prob_Q = symcfg.prob_Q
        kl[key].append( KLd(gibbs_P, prob_Q) )
        
        bv.cmp_graph(gibbs_P, prob_Q, ax=axis[i*2], label=key+str(j))
        printKL(arange(t), kl[key], ax=axis[i*2+1])


    plt.draw()
#    plt.show()
    fname = "_".join([symcfg.prefix for symcfg in symconfigs])
    plt.savefig('/home/shiogai/dir-png/%s_%d.png'%(fname,j))
    j+=1
    if not line:
        for symcfg in enumerate(symconfigs):
            symcfg.close()
        break




plt.show()
