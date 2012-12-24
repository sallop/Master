#!/usr/bin/python
#-*- coding: utf-8 -*-
from optparse import OptionParser


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


if __name__ == '__main__':
    print options
    print args
