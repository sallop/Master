#!/usr/bin/env python

from pylab import *

import Gnuplot

def sigma(x):
    return 1./(1.+exp(-x))

def prop(u,temperature):
    return sigma(u/temperature)

gp = Gnuplot.Gnuplot()

def fn(u,t):
    return "1./(1.+exp(-%(u)s/%(t)s))"%{"u":u,"t":t}
#    return "\frac{1}{1+\exp(-%(u)s/%(t)s)}"%{"u":u,"t":t}


gp(
"""
set xrange[-30:30]
set grid
set key left top
set terminal postscript enhanced color eps
set output 'prob_fix_u.eps'
set title '%(title)s'
f(u,t) = 1./(1.+exp(-u/t))
plot \
f(1,x) title '%(1)s',\
f(2,x) title '%(2)s',\
f(3,x) title '%(3)s',\
f(4,x) title '%(4)s',\
f(5,x) title '%(5)s',\
f(6,x) title '%(6)s',\
f(7,x) title '%(7)s',\
f(8,x) title '%(8)s',\
f(9,x) title '%(9)s'
"""%{"title":fn('u','T'),
     "1":fn(1,'T'),
     "2":fn(2,'T'),
     "3":fn(3,'T'),
     "4":fn(4,'T'),
     "5":fn(5,'T'),
     "6":fn(6,'T'),
     "7":fn(7,'T'),
     "8":fn(8,'T'),
     "9":fn(9,'T')}
)

gp("""
set xrange[-30:30]
set grid
set key left top
set terminal postscript enhanced color eps
set output 'prob_fix_T.eps'
set title '%(title)s'
f(u,x) = 1./(1.+exp(-u/x))
plot \
f(x,1) title '%(1)s',\
f(x,2) title '%(2)s',\
f(x,3) title '%(3)s',\
f(x,4) title '%(4)s',\
f(x,5) title '%(5)s',\
f(x,6) title '%(6)s',\
f(x,7) title '%(7)s',\
f(x,8) title '%(8)s',\
f(x,9) title '%(9)s'
"""%{"title":fn('u','T'),
     "1":fn('u',1),
     "2":fn('u',2),
     "3":fn('u',3),
     "4":fn('u',4),
     "5":fn('u',5),
     "6":fn('u',6),
     "7":fn('u',7),
     "8":fn('u',8),
     "9":fn('u',9),})
