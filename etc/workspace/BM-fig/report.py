#!/usr/bin/python
# -*- coding: utf-8 -*-

def figure_parallel2(fig1, title1, fig2, title2):
    print """\
\\begin{figure}[htbp]
 \\begin{center}
  \\begin{tabular}{cc}
   \\begin{minipage}{0.5\\hsize}
    \\includegraphics[scale=0.5]{%(fig1)s}
    \\caption{"%(title1)s"}
   \\end{minipage}
   &
   \\begin{minipage}{0.5\\hsize}
    \\includegraphics[scale=0.5]{%(fig2)s}
    \\caption{"%(title2)s"}
   \\end{minipage}
  \\end{tabular}
 \\end{center}
\\end{figure}"""%locals()

def figure_parallel3(fig1, title1, fig2, title2, fig3, title3):
    print """\
\\begin{figure}[htbp]
 \\begin{center}
  \\begin{tabular}{ccc}
   \\begin{minipage}{0.333\\hsize}
    \\includegraphics[scale=0.5]{%(fig1)s}
    \\caption{"%(title1)s"}
   \\end{minipage}
   &
   \\begin{minipage}{0.333\\hsize}
    \\includegraphics[scale=0.5]{%(fig2)s}
    \\caption{"%(title2)s"}
   \\end{minipage}
   &
   \\begin{minipage}{0.333\\hsize}
    \\includegraphics[scale=0.5]{%(fig3)s}
    \\caption{"%(title3)s"}
   \\end{minipage}
  \\end{tabular}
 \\end{center}
\\end{figure}"""%locals()



def simulation_condition(title, label, d, p, k, m, n, t, alpha):
    print """\
\\begin{table}
 \\begin{center}
  \\caption{"%(title)s"}
   \\begin{tabular}{|c|c|c|c|c|c|c|}
    \\hline
可視素子数 & 隱れ素子数 & k-gibbs & モデル標本数 & 訓練データ数 & 訓練時間 & 学習係数 \\\\
    \\hline
%(d)s    & %(p)s    & %(k)s   & %(m)s      & %(n)s     & %(t)s  & %(alpha)s\\\\
    \\hline
   \\end{tabular}
   \\label{tbl:%(label)s}
  \\end{center}
\\end{table}
"""%locals()

condition1 = {
    "title": "全結合Boltzmann Machine",
    "label": "FBM",
    "d"    : 3,
    "p"    : 1,
    "k"    : 5,
    "m"    : 1000,
    "n"    : 1000,
    "t"    : 1000,
    "alpha": 1.0
}

condition2 = {
    "title": "隱れ層有",
    "label": "GBM",
    "d"    : 3,
    "p"    : 1,
    "k"    : 5,
    "m"    : 1000,
    "n"    : 1000,
    "t"    : 1000,
    "alpha": 1.0
}

condition3 = {
    "title": "Restricted Boltzmann Machine",
    "label": "RBM",
    "d"    : 3,
    "p"    : 1,
    "k"    : 5,
    "m"    : 1000,
    "n"    : 1000,
    "t"    : 1000,
    "alpha": 1.0
}

simulation_figure1 = {
    "fig1"   : "dir-FBMBias/t-kl.dat",
    "fig2"   : "dir-FBMNoBias/t-kl.dat",
    "title1" : "全結合 バイアス有",
    "title2" : "全結合 バイアス無",
}

simulation_figure2 = {
    "fig1"   : "dir-RBMBias/t-kl.dat",
    "fig2"   : "dir-RBMNoBias/t-kl.dat",
    "title1" : "RBM バイアス有",
    "title2" : "RBM バイアス無",
}

simulation_GBM_Bias = {
    "fig1"   : "dir-GBMBiasIterate/t-kl.dat"  ,
    "fig2"   : "dir-GBMBiasNoIterate/t-kl.dat",
    "title1" : "GBM eq-",
    "title2" : "GBM eq-",
}

simulation_GBM_Iterate = {
    "fig1"   : "dir-GBMBiasIterate2/t-kl.dat",
    "fig2"   : "dir-GBMBiasIterate3/t-kl.dat",
    "title1" : "GBM $a_{j}$ + $\sum_{j} b_{j} h_{j}$",
    "title2" : "GBM $\sum_{i} a_{i} v_{i} \sum_{j} b_{j} h_{j}$",
}


simulation_figure5 = {
    "fig1"   : "dir-GBMBias/t-kl.dat",
    "fig2"   : "dir-GBMNoBias/t-kl.dat",
    "title1" : "GBM バイアス有",
    "title2" : "GBM バイアス無",
}

simulation_figure6 = {
    "fig1"   : "dir-GBMBias/t-kl.dat",
    "fig2"   : "dir-GBMNoBias/t-kl.dat",
    "title1" : "GBM バイアス有",
    "title2" : "GBM バイアス無",
}


#simulation_condition(**condition)
