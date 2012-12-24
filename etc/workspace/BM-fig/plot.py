#!/usr/bin/python
# -*- coding: utf-8 -*-

import subprocess
import os
import sys

def plot_command_fix(title, odir, file1, title1, file2, title2):
    return """\
set grid
set yrange[0:1.0]
set terminal postscript enhanced color eps
set title  "%(title)s" 
set output "%(odir)s/%(title)s.dat"
plot "%(file1)s" w l title "%(title1)s",\
     "%(file2)s" w l title "%(title2)s"
"""%locals()

def plot_command_var(title, odir, file1, title1, file2, title2):
    return """\
set grid
set terminal postscript enhanced color eps
set title  '%(title)s'
set output '%(odir)s/%(title)s.eps'
plot "%(file1)s" w l title "%(title1)s",\
     "%(file2)s" w l title "%(title2)s"
"""%locals()

def plot_command_var4(title, odir,
                      file1, title1,
                      file2, title2,
                      file3, title3,
                      file4, title4):
    return """\
set grid
set terminal postscript enhanced color eps
set title  '%(title)s'
set output '%(odir)s/%(title)s.eps'
plot "%(file1)s" w l title "%(title1)s",\
     "%(file2)s" w l title "%(title2)s",\
     "%(file3)s" w l title "%(title3)s",\
     "%(file4)s" w l title "%(title4)s"
"""%locals()


simulation_FBM = {
    "title"  : "FBM-cmp-bias",
    "odir"   : "dir-fig",
    "file1"   : "dir-FBMBias/t-kl.dat",
    "file2"   : "dir-FBMNoBias/t-kl.dat",
    "title1" : "FBM bias  on",
    "title2" : "FBM bias off",
}

simulation_RBM = {
    "title"  : "RBM-cmp-bias",
    "odir"   : "dir-fig",
    "file1"   : "dir-RBMBias/t-kl.dat",
    "file2"   : "dir-RBMNoBias/t-kl.dat",
    "title1" : "RBM bias  on",
    "title2" : "RBM bias off",
}

simulation_GBM_Iterate = {
    "title"  : "GBM-cmp-iterate",
    "odir"   : "dir-fig",
    "file1"   : "dir-GBMNoBiasIterate/t-kl.dat"  ,
    "file2"   : "dir-GBMNoBiasNoIterate/t-kl.dat",
    "title1" : "GBM    iterate",
    "title2" : "GBM no iterate",
}

simulation_GBM_Bias = {
    "title"  : "GBM-cmp-bias",
    "odir"   : "dir-fig",
    "file1"   : "dir-GBMBiasIterate2/t-kl.dat",
    "file2"   : "dir-GBMBiasIterate3/t-kl.dat",
    "title1" : "GBM bias part a_i b_j",
    "title2" : "GBM bias all  a_i b_j",
}

simulation_GBM_Iterate = {
    "title"  : "GBM-cmp-iterate",
    "odir"   : "dir-fig",
    "file1"   : "dir-GBMNoBiasIterate/t-kl.dat"  ,
    "file2"   : "dir-GBMNoBiasNoIterate/t-kl.dat",
    "title1" : "GBM    iterate",
    "title2" : "GBM no iterate",
}

simulation_GBM_all = {
    "title" : "GBM-cmp-all",
    "odir"  : "dir-fig",
    "file1" : "dir-GBMNoBiasIterate/t-kl.dat"  ,
    "file2" : "dir-GBMNoBiasNoIterate/t-kl.dat",
    "file3" : "dir-GBMBiasIterate2/t-kl.dat",
    "file4" : "dir-GBMBiasIterate3/t-kl.dat",
    "title1": "GBM           iterate",
    "title2": "GBM        no iterate",
    "title3": "GBM bias part a_i b_j",
    "title4": "GBM bias all  a_i b_j"
}

simulation_BM_bias = {
    "title" : "BM-bias",
    "odir"  : "dir-fig",
    "file1" : "dir-FBMBias/t-kl.dat"  ,
    "title1": "FBM",
    "file2" : "dir-RBMBias/t-kl.dat",
    "title2": "RBM",
    "file3" : "dir-GBMBiasIterate2/t-kl.dat",
    "title3": "GBM bias part a_i b_j",
    "file4" : "dir-GBMBiasIterate3/t-kl.dat",
    "title4": "GBM bias all  a_i b_j"
}

simulation_BM_nobias = {
    "title"  : "BM-nobias",
    "odir"   : "dir-fig",
 
    "file1"  : "dir-FBMNoBias/t-kl.dat",
    "title1" : "FBM",
    
    "file2"  : "dir-RBMNoBias/t-kl.dat",
    "title2" : "RBM",

    "file3"  : "dir-GBMNoBiasIterate/t-kl.dat",
    "title3" : "GBM    iterate",

    "file4"  : "dir-GBMNoBiasNoIterate/t-kl.dat",
    "title4" : "GBM no iterate"
}


plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var(**simulation_FBM))
plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var(**simulation_RBM))
plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var(**simulation_GBM_Bias))
plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var(**simulation_GBM_Iterate))

plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var4(**simulation_GBM_all))

plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var4(**simulation_BM_bias))

plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
plot.communicate(plot_command_var4(**simulation_BM_nobias))
