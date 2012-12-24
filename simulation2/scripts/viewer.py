#!/usr/bin/ipython
#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
import pylab
import time
from math import *
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator

def print_KL(time,kl,ax=None):
    """
    print KL divergence
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    ax.cla()
    ax.set_xlim(   0,  100)
    ax.set_ylim(-0.1,  1.1)
    ax.plot(time,kl)
    #print "ax.ylim()=",ax.get_ylim()

def cmp_graph(P,Q,ax=None,label=None):
    """
    compare distribution
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    ax.cla()
    ind   = np.arange(len(P))
    width = 1.0/2.0 - 0.1
    rects1 = ax.bar(ind      , P, width, color='r')
    rects2 = ax.bar(ind+width, Q, width, color='y')
    #ax.set_ylabel('probability')
    ax.set_title (label)

    ax.set_xticks(ind+width)
    ax.set_xticklabels(tuple(["%s"%(i) for i,p in enumerate(P)]))
    #ax.legend((rects1[0],rects2[0]),('P','Q'))
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()

    autolabel(rects1)
    autolabel(rects2)
    ax.set_ylim(0.0, 1.0)
    ax.set_label(label)



def cmp_energy(ene, ax):
    """
    compare energy
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    n = len(ene)
    ind = np.arange(n)
    width = 0.35
    ax.cla()
    ax.bar(ind, ene, width, color = 'r')
    ax.set_xlim(0,n)
    ax.set_ylabel('energy')
    return

def hinton(L, maxWeight=None, ax=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if not maxWeight:
        maxWeight = 2**np.ceil(np.log(np.abs(L).max())/np.log(2))

    #ax.cla()
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal','box')
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_locator(NullLocator())

    for (x,y),w in np.ndenumerate(L):
        if w > 0 : color = 'white'
        else     : color = 'black'
        size = np.sqrt(abs(w))
        rect = Rectangle([x-size/2,y-size/2],
                         size,
                         size,
                         facecolor=color,
                         edgecolor=color)
        ax.add_patch(rect)
    ax.autoscale_view()
    # Reverse the yaxis limits
    ax.set_ylim(*ax.get_ylim()[::-1])

def torus_plot(xc,yc,zc,rd,n,ax=None):
    if ax == None:
        fig = plt.figure()
        rect = fig.add_subplot(111).get_position()
        ax = Axes3D(fig,rect)

    zc2 = zc + rd
    # Visible Layer Circle
    xs1,ys1,zs1 = Circle3D(xc,yc,zc,rd).create(n)
    xs1 = np.array(xs1 + xs1[:1])
    ys1 = np.array(ys1 + ys1[:1])
    zs1 = np.array(zs1 + zs1[:1])
    ax.plot(xs1,ys1,zs1)
    # Hidden Layer Circle
    xs2,ys2,zs2 = Circle3D(xc,yc,zc2,rd).create(n)
    xs2 = np.array(xs2 + xs2[:1])
    ys2 = np.array(ys2 + ys2[:1])
    zs2 = np.array(zs2 + zs2[:1])
    ax.plot(xs2,ys2,zs2)
    # Visible to Hiden Layer Bar
    xs3 = [(x1, x2) for x1,x2 in zip(xs1,xs2)]
    ys3 = [(y1, y2) for y1,y2 in zip(ys1,ys2)]
    zs3 = [(zc,zc2) for i in xrange(len(xs3))]
    for x3,y3,z3 in zip(xs3,ys3,zs3):
        ax.plot(x3,y3,z3)
