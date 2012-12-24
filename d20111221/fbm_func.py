#!/usr/bin/env python
#-*- coding: utf-8 -*-
from pylab import *

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
        for j in xrange(i,len(L[0])):
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
    f = lambda p,q: p*log(p) - p*log(q)
    for p,q in zip(P,Q): d += f(p,q)
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

