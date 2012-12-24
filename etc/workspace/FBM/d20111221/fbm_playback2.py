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
axis = [fig.add_subplot(211),
        fig.add_subplot(212),]

kl_result = []

for v in all_v:
    print "v = ", v

for line in open(files['t_gibbs'],"r"):
    import time
    t, gibbs_P = read_t_gibbs(line)
    kl = KLd(gibbs_P, prob_Q)
    kl_result.append(kl)

    print "t = ",t,
    for i,p in enumerate(gibbs_P):
        print "p(x=%d)=%4.3f"%(i,p),

    print "kl=%4.3f"%(kl)
    bv.cmp_graph(gibbs_P,prob_Q,ax=axis[0])
    bv.print_KL  (arange(t), kl_result, ax = axis[1])

    plt.draw()

plt.show()
