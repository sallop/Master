#!/usr/bin/env python
#-*- coding: utf-8 -*-
from pylab import *
import viewer as bv
import matplotlib.pyplot as plt
import hinton_dia1

from fbm_func import *
from fbm_option import *
from readcfg import *
from read_ublas import *

dir = "dir-dat"
prefix = sys.argv[1] if sys.argv[1] else "K23x3"
identifies = ["cfg","mat","t_gibbs"]
suffix = ".dat"
files = dict((id, "%s/%s-%s%s"%(dir,prefix,id,suffix))
             for id in identifies)

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


symcfg = {}
fr = open(files['cfg'],'r')
for i,line in enumerate(fr):
    if i > 32: break
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
g_annealing_schedule_idx = int(symcfg['g_annealing_schedule_idx'].pop())
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
print "g_annealing_schedule_idx=",g_annealing_schedule_idx
print "alpha=",alpha
print "prob_Q=",prob_Q
print "all_v=",all_v 
#assert False

plt.ion()
fig = plt.figure()
axis = [fig.add_subplot(411),
        fig.add_subplot(412),
        fig.add_subplot(413),
        fig.add_subplot(414),]
axis[2].set_xlim(-1,16)
kl_result = []

for line in open(files['mat'],"r"):
    import time
    t, mat = read_t_mat(line)
    L = array(mat, dtype=float)
    prob_P = [boltzmann_distribution(v,L,all_v) for v in all_v]
    ene_result = map(lambda v : energy(v,L), all_v)
    kl = KLd(prob_P, prob_Q)
    kl_result.append(kl)

    for i,v in enumerate(all_v):
        print "v[%d]="%(i),v

    for i,ene in enumerate(ene_result):
        print "E%d=%8.4f"%(i,ene),
    print "kl=%8.4f"%(kl)
    hinton_sub   (L, ax = axis[0])
    bv.cmp_energy(ene_result, ax = axis[1])
    bv.print_KL  (arange(t), kl_result, ax = axis[2])
    bv.cmp_graph(prob_P, prob_Q, ax = axis[3])
    plt.draw()

plt.show()
