#!/usr/bin/python
#-*- coding: utf-8 -*-

from math import *

class drawing_object:
    pass

class Circle2D(drawing_object):
    """
    drawing circle object for edges
    """
    def __init__(self,x,y,r):
        self.xc = x
        self.yc = y
        self.rd = r
        
    def get(self,xs,ys):
        rd = self.rd
        n = len(xs)
        dt = 2.*pi/n
        for i in xrange(n):
            xs[i] = rd*cos(i*dt)
            ys[i] = rd*sin(i*dt)
        return xs, ys

    def create(self,n):
        rd = self.rd
        dt = 2.*pi/n
        xs = [rd*cos(i*dt) for i in xrange(n)]
        ys = [rd*sin(i*dt) for i in xrange(n)]
        return (xs, ys)

class Circle3D(drawing_object):
    """
    drawing circle object for edges at 3D
    """
    def __init__(self,x,y,z,r):
        self.xc = x
        self.yc = y
        self.zc = z
        self.rd = r

    def set_core(self,x,y,z):
        self.xc = x
        self.yc = y
        self.zc = z

    def set_radius(self,r):
        self.rd = r

    def get_core(self):
        return (self.xc, self.yc, self.zc)

    def get(self,xs,ys,zs):
        rd = self.rd
        zc = self.zc
        n  = len(xs)
        dt = 2.*pi/n
	zs = [zc for i in xrange(n)]
        for i in xrange(n):
            xs[i] = rd*cos(i*dt)
            ys[i] = rd*sin(i*dt)

        return xs, ys, zs

    def create(self,n):
        dt = 2.*pi/n
        rd = self.rd
        zc = self.zc
        xs = [rd*cos(i*dt) for i in xrange(n)]
        ys = [rd*sin(i*dt) for i in xrange(n)]
        zs = [zc           for i in xrange(n)]
        return xs, ys, zs

if __name__ == '__main__':
    import Gnuplot
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.axes3d import Axes3D

    circle = Circle2D(2,3,5)
    g = Gnuplot.Gnuplot()
    xs, ys ,zs = circle.create(12)
    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    xs2 = [0 for i in xrange(12)]
    ys2 = [0 for i in xrange(12)]
    circle.get(xs2,ys2)

    fig = plt.figure()
    # figure 1
    rect= fig.add_subplot(2,1,1).get_position()
    ax = Axes3D(fig, rect)
    ax.scatter(xs,ys,zs)
    # figure 2
    rect= fig.add_subplot(2,1,2).get_position()
    ax = Axes3D(fig, rect)
    ax.plot(xs,ys,zs)


    plt.show()
#     g("plot '-' with points, '-' with points")
#     for x,y in zip(xs,ys):
#         g('%f %f'%(x,y))
#     g("e")
#     for x,y in zip(xs,ys):
#         g('%f %f'%(x,y))
#     g("e")
#     g("pause 2.5")
