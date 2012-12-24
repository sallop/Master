#!/usr/bin/ipython
#-*- coding: utf-8 -*-
from pylab import *
import time
from math import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator
from viewer import *

def dtob(n, size):
    ret = np.zeros(size)
    for i in xrange(ret.size):
        ret[size-i-1]=(n>>i)&0x01
    return ret

def btod(x):
    value=0.0
    for i in xrange(x.size):
        value += x[x.size-i-1]*(0x01 << i)
    return value

def sigma(u): return 1./(1.+exp(-u))

def energy(v, h, W, a, b):
    term1, term2, term3 = 0.0, 0.0, 0.0
    D = v.size
    P = h.size
    term1 = dot(    a, v)
    term2 = dot(    b, h)
    term3 = dot(  v.T, W)
    term3 = dot(term3, h)
    return - term1 - term2 - term3

def propup(j, W, v, b):
    return sigma(b[j] + dot(W[:,j],v))

def propdown(i, W, h, a):
    return sigma(a[i] + dot(W[i,:],h))

def partition_function(W,a,b,all_v,all_h):
    ret = 0.0
    for v in all_v:
        for h in all_h:
            ret += exp(-energy(v,h,W,a,b))
    return ret

def mean_field_update(j,W,v):
    return sigma(dot(v,W[:,j]))
# sigma(b[j]+np.dot(v,W[:,j]))

def KLd(P,Q):
    d = 0.0
    for p,q in zip(P,Q) : d += p*log(p/q)
    return d

def prob_model_assign_v(v,W,a,b,all_v,all_h):
    numer, denom = 0.0, 0.0
    for h in all_h : numer += exp(-energy(v,h,W,a,b))
    denom = partition_function(W,a,b,all_v,all_h)
    return numer/denom

D = 3
P = 2
K = 50
M = 1000
N = 1000
T = 1000
alpha = 1.0                     # learning rate

probP = zeros(0x01<<D)
probQ = array([0.10, 0.10, 0.05, 0.05, 0.10, 0.10, 0.40, 0.10])
accuQ = add.accumulate(probQ)
accuQ[-1]=1.0

all_v = [dtob(n,D) for n in xrange(0x01<<D)]
all_h = [dtob(n,P) for n in xrange(0x01<<P)]

v = rand(D)
a = rand(D)
vf= rand(D)
h = rand(P)
hf= rand(P)
b = rand(P)
mu= rand(P)
W       = rand(D,P)
W_inc   = rand(D,P)
W_data  = rand(D,P)
W_model = rand(D,P)

kl = 0.0

plt.ion()
fig = plt.figure()

for t in xrange(1,T):
    # get a DATA dependent expection
    W_data[:]=0
    for n in xrange(N):
        # find a small value
        rnd = rand()
        for k,aq in enumerate(accuQ):
            if rnd < aq: break

        vk = all_v[k]           # not side effect

        mu = rand(P)
        for i in xrange(100):
            j = randint(0,P-1)
            mu[j] = mean_field_update(j,W,vk)

        W_data += outer(vk,mu)

    # get MODEL dependent expection
    W_model[:]=0
    for m in xrange(M):
        for k in xrange(K):
            j = randint(0,P-1)
            i = randint(0,D-1)
            rnd = rand()
            hf[j] = 1 if propup  (j,W,v,b) > rnd else 0
            vf[i] = 1 if propdown(i,W,h,a) > rnd else 0
        W_model += outer(vf,hf)
    
    W_data  /= N
    W_model /= M
    W += alpha*(W_data - W_model) # W is not symmetric matrix
    #symmetrize_U(W)
#    alpha = 1.0/t # descrese learning rate

    # Print KL divergence
    for i,v in enumerate(all_v):
        probP[i] = prob_model_assign_v(v,W,a,b,all_v,all_h)

    print "kl=%8.4f\t"%(KLd(probP,probQ)),
    for p in probP:
        print "%8.4f "%(p),
    print
    plt.clf()
    axis = [fig.add_subplot(1,3,1),
            fig.add_subplot(1,3,2),
            fig.add_subplot(1,3,3),]
    axis[0].set_title("W")
    axis[1].set_title("W data dependent")
    axis[2].set_title("W model dependent")
    hinton(W      , ax=axis[0])
    hinton(W_data , ax=axis[1])
    hinton(W_model, ax=axis[2])
    plt.draw()


