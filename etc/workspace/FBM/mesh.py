#!/usr/bin/env python
#-*- coding: utf-8 -*-
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import Gnuplot
#import matplotlib.pyplot as plt

gp = Gnuplot.Gnuplot()

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


def rotate_matrix(x,phi=0,theta=0,psi=0):
    R1 = array([[cos(phi), -sin(phi), 0],
                [sin(phi),  cos(phi), 0],
                [       0,         0, 1]])

    R2 = array([[1,          0,           0],
                [0, cos(theta), -sin(theta)],
                [0, sin(theta),  cos(theta)]])
    R3 = array([[cos(psi), -sin(psi), 0],
                [sin(psi),  cos(psi), 0],
                [       0,         0, 1]])

    return dot(R3,dot(R2,dot(R1,x)))

def polar(r,theta,phi):
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    return x,y,z

def vertical_circle(num,phi,r=1.0):
    f = lambda theta: polar(r,theta,phi)
    return map(f, linspace(0., 2.*pi, num+1))

def horizon_circle(num,theta,r=1.0):
    f = lambda phi: polar(r, theta, phi)
    return map(f, linspace(0., 2.*pi, num+1))


def assert_circle_type(c):
    print "type(c)=",type(c)
    print "type(c[0])=",type(c[0])
    print "type(c[0][0])=",type(c[0][0])
    print c
    for xyz in c:
        print xyz

    print c[0][0]


def mesh(M):
    r = 1.0
    N = len(M)
    width = ceil(sqrt(D))
    height= N/width
    xs = linspace(0, width  , width  )
    ys = linspace(0, N/width, N/width)
    pos = [(x,y) for x in xs for y in ys]
    gp("""
set grid
set xrange[-0.5:4.5]
set yrange[-0.5:4.5]
""")

    gp("plot '-' with linespoint 1, '-' with points 3")
    for src in xrange(N):
        for dst in xrange(N):        
            if M[src,dst] > 0.0:
                gp("%f %f"%pos[src])
                gp("%f %f"%pos[dst])

            gp("\n")
    gp("e")

    for idx in xrange(N):
        if M[idx,idx] >= 2.0:
            gp("%f %f"%pos[idx])
            gp("%f %f"%pos[idx])

        gp("\n")
    gp("e")


    for idx in xrange(N):
        print "idx=%d\t"%(idx),pos[idx]

    gp("pause 10.0")
    return


def doughnut(vcircles, a1, a2):
    hcircles = zip(*vcircles)
    gp("""
set view %(a1)f, %(a2)f
set xrange[-2.:2.]
set yrange[-2.:2.]
set zrange[-2.:2.]
"""%{"a1":a1,"a2":a2})

    gp("splot '-' with lines")
    for vc in vcircles:
        for p in vc:
            gp("%f\t%f\t%f"%tuple(p))
        gp("\n")

    for hc in hcircles:
        for p in hc:
            gp("%f\t%f\t%f"%tuple(p))
        gp("\n")
    
    gp("e")
    gp("pause 0.1")
    return


def doughnut2(W, a1, a2):
    gp("""
set view %(a1)f, %(a2)f
set grid
set xrange[-2.:2.]
set yrange[-2.:2.]
set zrange[-2.:2.]
"""%{"a1":a1,
     "a2":a2})

    # vertical circle iterator index i, max D
    # horizon circle iterator index j, max P
    D = len(W)
    P = len(W[0])
    vr,hr = 0.5,1.5
    
    hc_offset = horizon_circle(D, 0.5*pi, hr)
    phis = linspace(0.,2.*pi,D+1)

    vcircles = []
    for p0,phi in zip(hc_offset,phis):
        vc = vertical_circle(P,phi,vr)
        vc = vc + array(p0)
        vcircles.append(vc)

    hcircles = zip(*vcircles)
            
    gp("splot '-' with lines")
    for i,vc in enumerate(vcircles):
        for j in xrange(P):
            if W[i][j] > 0:
                gp("%f\t%f\t%f\n"%tuple(vc))
            elif W[i][j] <= 0:
                gp("\n")
                
    for i,hc in enumerate(hcircles):
        for j in xrange(P):
            if W[i][j] > 0:
                gp("%f\t%f\t%f\n"%tuple(hc))
            elif W[i][j] <= 0:
                gp("\n")

    gp("e")
    gp("pause 3.2")
    return

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

def make_kneighbor_onbit(xs,k):
    dsize = len(xs[0])
    for i,x in enumerate(xs):
        x[:]=0
        for j in xrange(k):
            x[(i+j)%dsize]=1
    return xs

D = 16
P = 16
L = zeros(shape=(D,D),dtype=float)
W = zeros(shape=(D,P),dtype=float)
J = zeros(shape=(P,P),dtype=float)
vr = 0.5
hr = 1.5

all_v = [zeros(D) for i in xrange(D)]
make_kneighbor_onbit(all_v,2)

#L = reduce(lambda a,b:a+outer(b,b),all_v)
for v in all_v:
    L += outer(v,v)
print L

mesh(L)


phi, theta = 0.0, 0.0

points = horizon_circle(D,0.5*pi, hr)
phis = linspace(0.,2.*pi,D+1)
# vcircles = [vertical_circle(P,phi,0.3) + array(p0)
#             for p0,phi in zip(points,phis)]
vcircles = []
for p0,phi in zip(points,phis):
    vc = vertical_circle(P,phi, vr)
    vc = vc + array(p0)
#    assert_circle_type(vc)
    vcircles.append(vc)



# 垂直方向
# for i,vc in enumerate(vcircles):
#     print type(vc)
#     print vc

#     print type(vc[0])
#     print vc[0]

#     print type(vc[0][0])
#     print vc[0][0]

# #    xs,ys,zs = zip(*vc)
#     for p in vc:
#         gp("%f\t%f\t%f"%tuple(p))

#     gp("\n")


# 水平方向
# print type(vcircles)
# print vcircles
# print type(vcircles[0])
# print vcircles[0]
# print type(foo)
# print foo

# hcircles = zip(*vcircles)
# for i,hc in enumerate(hcircles):
#     for p in hc:
#         gp("%f\t%f\t%f"%tuple(p))

#     for z in zs:
#         for y in ys:
#             for x in xs:
#                 gp("%f\t%f\t%f"%(x,y,z))
#                 print "%f\t%f\t%f"%(x,y,z)
#         gp("\n")
#         print "\n"

# gp("e")
# gp("pause 3.2")
        
# for a1 in linspace(0,180,D+1):
#     for a2 in linspace(0,180,P+1):
#         #doughnut2(L, a1, a2)
#         doughnut(vcircles, a1, a2)



# for p in points: gp("%f\t%f\t%f"%p)
# print "draw vertical=",len(phis)
# for vc in vcircles:
#     print "vc"
#     for p in vc:
#         gp("%f\t%f\t%f\n"%tuple(p))
#         print "%f\t%f\t%f"%tuple(p)
#     gp("\n")
# gp("e")
# gp("pause 10.0")
# fig = figure()
# rect = fig.add_subplot(111).get_position()
# ax = Axes3D(fig,rect)
# ax.set_xlim((-2,2))
# ax.set_ylim((-2,2))
# ax.set_zlim3d((-2,2))
# xs,ys,zs=zip(*points)
# ax.plot(xs,ys,zs)
# for vc in vcircles:
#     xs,ys,zs = zip(*vc)
#     ax.plot(xs,ys,zs)
# show()
