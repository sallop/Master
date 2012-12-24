#!/usr/bin/env python
#-*- coding: utf-8 -*-
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import Gnuplot
from draw_object import *
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
    

def torus_gp(**keys):
    gp("""
set grid
set view %(a1)f, %(a2)f
set xrange[-2:2]
set yrange[-2:2]
set terminal postscript enhanced color eps
set output "dir-pudding/%(output)s.eps"
"""%keys)
       #{"a1":a1,"a2":a2})
    D = 16
    P = 8
    hr1 = 1.5
    hc1 = zip(*Circle3D(0., 0., 0., hr1).create(D))
    hr2 = 1.0
    hc2 = zip(*Circle3D(0., 0., 10., hr2).create(P))
    
    phis1 = linspace(0.,2.*pi, D)
    phis2 = linspace(0.,2.*pi, P)
    
    gp("splot '-' with lines")
    for p1 in hc1:
        gp("%f %f %f"%(p1))
        pass
    gp("\n")

    for p2 in hc2:
        gp("%f %f %f"%(p2))
        pass
    gp("\n")

    for p1 in hc1:
        for p2 in hc2:
            gp("%f %f %f"%(p1))
            gp("%f %f %f"%(p2))
            pass
        pass
    gp("e")
    return

def mesh(W):
    r = 1.0
    D,P = len(W),len(W[0])
    N = D+P
    xs = linspace(0, P, P)
    ys = linspace(0, D, D)
#     x = r*cos(linspace(0, 2.0*pi, N))
#     y = r*sin(linspace(0, 2.0*pi, N))
    npoints = []
    for x in xs:
        for y in ys:
            npoints.append((x,y))

    gp("""
set grid
#set terminal postscript enhanced eps color
#set output 'dir-mesh/%(output)s'
""")

    gp("plot '-' with lines")
    for n in xrange(N):
        i = n/D
        j = n%D
        if W[i,j] > 0.0:
            print W[i,j]
            gp("%d %d"%(xs[i],ys[j]))
        else:
            gp("\n")

    gp("\n")

    for n in xrange(N):
        i = n%D
        j = n/D
        if W[i,j] > 0.0:
            print W[i,j]
            gp("%d %d"%(xs[i],ys[j]))
        else:
            gp("\n")

    gp("\n")
    gp("e")

    gp("pause 30.0")
    return


def doughnut(vcircles, a1, a2):
    hcircles = zip(*vcircles)
    gp("""
set view %(a1)f, %(a2)f
set grid
#set terminal postscript enhance color eps
#set output "doughnut.eps"
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


def doughnut2(W, **keys):
    gp("""
set view %(a1)f, %(a2)f
set grid
set xrange[-2.:2.]
set yrange[-2.:2.]
set zrange[-2.:2.]
set terminal postscript enhanced color eps
set output 'dir-doughnut/%(output)s.eps'
"""%keys)

    # vertical circle iterator index i, max D
    # horizon circle iterator index j, max P
    D = len(W)
    P = len(W[0])
    vr,hr = 0.5,1.5
    
    hc_offset = horizon_circle(D, 0.5*pi, hr)
    phis = linspace(0.,2.*pi,D)

    vcircles = []
    for p0,phi in zip(hc_offset,phis):
        vc = vertical_circle(P,phi,vr)
        vc = vc + array(p0)
        vcircles.append(vc)

    hcircles = zip(*vcircles)
            
    gp("splot '-' with lines")
    for i,vc in enumerate(vcircles):
        for j,p in enumerate(vc):
            gp("%f\t%f\t%f"%tuple(p))
        gp("\n")

    for i,hc in enumerate(hcircles):
        for j,p in enumerate(hc):
            gp("%f\t%f\t%f"%tuple(p))


    gp("e")


def D2toD1(W, **keys):
    gp("""
set view %(a1)f, %(a2)f
set grid
set xrange[-3.:3.]
set yrange[-3.:3.]
set zrange[-3.:3.]
set terminal postscript enhanced color eps
set output 'dir-D2toD1/%(output)s.eps'
"""%keys)

    # vertical circle iterator index i, max D
    # horizon circle iterator index j, max P
    D = len(W)
    P = len(W[0])
    vr,hr,hr2 = 0.2, 1.5, 0.8
    
    hc_offset = horizon_circle(D, 0.5*pi, hr)

    phis = linspace(0.,2.*pi,D)

    vcircles = []
    for p0,phi in zip(hc_offset,phis):
        vc = vertical_circle(P,phi,vr)
        vc = vc + array(p0)
        vcircles.append(vc)

    hcircles = zip(*vcircles)
            
    gp("splot '-' with lines")
    for i,vc in enumerate(vcircles):
        for j,p in enumerate(vc):
            gp("%f\t%f\t%f"%tuple(p))
        gp("\n")

    for i,hc in enumerate(hcircles):
        for j,p in enumerate(hc):
            gp("%f\t%f\t%f"%tuple(p))

    hcircle2 = horizon_circle(D, 0.5*pi, hr2)
    hcircle2 = [p1+array((0.,0.,2.5)) for p1 in hcircle2]
    
    for i,p in enumerate(hcircle2):
        gp("%f\t%f\t%f"%tuple(p))

    gp("\n")
    for p1,p2 in zip(hc_offset, hcircle2):
        gp("%f\t%f\t%f"%tuple(p1))
        gp("%f\t%f\t%f"%tuple(p2))
        gp("\n")
    gp("e")



    # for i,vc in enumerate(vcircles):
#         for j in xrange(P):
#             if W[i][j] > 0:
#                 gp("%f\t%f\t%f\n"%tuple(vc[j]))
#             elif W[i][j] <= 0:
#                 gp("\n")
                
#     for i,hc in enumerate(hcircles):
#         for j in xrange(P):
#             if W[i][j] > 0:
#                 gp("%f\t%f\t%f\n"%tuple(hc[j]))
#             elif W[i][j] <= 0:
#                 gp("\n")


#    gp("pause 0.2")
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

#assert False

phi, theta = 0.0, 0.0

points = horizon_circle(D,0.5*pi, hr)
phis = linspace(0.,2.*pi,D+1)
vcircles = [vertical_circle(P,phi,0.3) + array(p0)
            for p0,phi in zip(points,phis)]
# vcircles = []
# for p0,phi in zip(points,phis):
#     vc = vertical_circle(P,phi, vr)
#     vc = vc + array(p0)
# #    assert_circle_type(vc)
#    vcircles.append(vc)



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

for i,a1 in enumerate(linspace(0,180,D+1)):
    for j,a2 in enumerate(linspace(0,180,P+1)):
        #torus_gp(a1=a1,a2=a2,output="%s"%(i+j))
        #doughnut2(L, a1=a1, a2=a2, output="%s"%(i+j))
        D2toD1(L, a1=a1, a2=a2, output="%s"%(i+j))
        #doughnut(vcircles, a1, a2)

#doughnut(vcircles, 60, 30)
gp("pause 3.2")

#mesh(L)

# for a1 in linspace(0, 180, D+1):
#     for a2 in linspace(0, 180, P+1):

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
