

def figure_parallel(fig1, title1, fig2, title2)
  temp=<<-EOS
\\begin{figure}[htbp]
 \\begin{center}
  \\begin{tabular}{cc}
   \\begin{minipage}{0.5\\hsize}
    \\includegraphics[scale=0.5]{#{fig1}}
    \\caption{#{title1}}
   \\end{minipage}
   &
   \\begin{minipage}{0.5\\hsize}
    \\includegraphics[scale=0.5]{#{fig2}}
    \\caption{#{title2}}
   \\end{minipage}
  \\end{tabular}
 \\end{center}
\\end{figure}
EOS
  return temp
end

def simulation_condition(title, label, d, p, k, m, n, t, alpha)
  temp=<<-EOS
\\begin{table}
 \\begin{center}
  \\caption{"#{title}"}
  \\begin{tabular}{|c|c|c|c|c|c|c|}
   \\hline
   可視素子数 & 隱れ素子 & k-gibbs & モデルの標本数 & 訓練データ & 訓練時間 & 学習率\\
   \\hline
   #{d} & #{p} & #{k} & #{m} & #{n} & #{t} & #{alpha}\\
   \\hline
  \\end{tabular}
  \\label{tbl:#{label}}
 \\end{center}
\\end{table}
EOS
end

puts figure_parallel('"dir-FBMBias/t-kl.dat"'  , '"FBM-バイアス有"',
                     '"dir-FBMNoBias/t-kl.dat"', '"FBM-バイアス無"')


puts figure_parallel('"dir-FBMBias/t-kl.dat"'  , '"FBM-バイアス有"',
                     '"dir-FBMNoBias/t-kl.dat"', '"FBM-バイアス無"')


puts simulation_condition("シミュレーション条件","sim-cnd",
                          3, 1, 5, 1000, 1000, 1000, 1)
