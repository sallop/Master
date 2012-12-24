#!/usr/bin/python
#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
from pylab import *
import time
import os
import sys
import subprocess
from math import *
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator
from optparse import OptionParser
import Gnuplot

parser = OptionParser()

parser.add_option('-m',
                  '--matrix',
                  dest   ='matrix',
                  action ='store_true',
                  default= False,
                  help   ='show a hinton diagram')

parser.add_option('-k',
                  '--kld',
                  action ='store_true',
                  dest   ='kld',
                  default= False,
                  help   ='don\'t print status messages to stdout')

parser.add_option('-t',
                  '--time',
                  action ='store',
                  nargs  = 2,
                  type   ='int',
                  dest   ='time',
                  default= False,
                  help   ='don\'t print status messages to stdout')


parser.add_option('-c',
                  '--compare',
                  action ='store_true',
                  dest   ='compare',
                  default= False,
                  help   ='compare probabilities')


parser.add_option('-b',
                  '--bar',
                  action ='store_true',
                  dest   ='bar',
                  default= False,
                  help   ='print probability histgram')

parser.add_option('-o',
                  '--output',
                  action ='store',
                  dest   ='output',
                  type   = 'string',
                  default= False,
                  help   ='output file name')

(options, args) = parser.parse_args()

def get_line (f, row):
    datas = get_lines(f,row,row+1)
    return datas

def get_lines(f, start, end):
    lines = f.readlines()
    datas = [map(lambda x:float(x),line.split()) for line in lines]
    f.seek(0,os.SEEK_SET)
    return datas[start:end]

def main(options, args):
    times = getattr(options,'time')
    
    for filename in args:
        fr = open(filename,'r')
        data = get_lines(fr,*times)
        dataT=zip(*data)
        t = dataT[0]
        kl= dataT[1]
        prob = dataT[2:]
        mats = zip(*[mat for mat in prob])
        # plot matrix
        if getattr(options,'matrix'):
            gnuplot('matrix',
                    *mats,
                    title='3x2',
                    rowsz=3,
                    colsz=2)
        #
        if getattr(options,'kld'):
            gnuplot('lines', t, kl, title='kld')

        if getattr(options,'compare'):
            gnuplot('lines', t, kl, title='compare')

        if getattr(options,'bar'):
            gnuplot('boxes', *prob, title='')

        fr.close()



if __name__ == '__main__':
    main(options, args)

