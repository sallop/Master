#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import sys
import re

re_tokens = [re.compile(k) for k in r"""
(D)=(\d+)
(P)=(\d+)
(K)=(\d+)
(M)=(\d+)
(N)=(\d+)
(T)=(\d+)
(alpha)=(.+)
(all_v\[i\])=(.+)
(prob_Q)=(.+)
(g_energy_idx)=(\d+)
(g_denergy_idx)=(\d+)
(g_prop_idx)=(\d+)
(g_annealing_schedule_idx)=(\d+)
""".split()]

#(all_v\[i\])=(.+)

def read_config(line, symcfg={}):
    for re_token in re_tokens:
        result = re_token.search(line)
        if result:
            (k,v) = result.groups()
            (s,e) = result.span()
            v = v.strip()
            if symcfg.has_key(k):
                symcfg[k].append(v)
            else:
                symcfg[k] = [v]
    return symcfg

if __name__ == '__main__':
    DIR = "dir-dat"
    prefix = "L3x3"
    identifies = ["cfg","mat","t_gibbs"]
    suffix = ".dat"
    files = dict((id,DIR +"/"+ prefix +"-"+ id + suffix) for id in identifies)

    symcfg = {}
    fr = open(files['cfg'],'r')
    for i,line in enumerate(fr):
        if i > 19: break
        read_config(line,symcfg)
    fr.close()

    for k,v in symcfg.iteritems():
        print "k=",k,
        print "v=",v
