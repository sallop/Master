#!/usr/bin/ruby

dirnames = ["../RBM/dir-dat","../GBM/dir-P2"]

IO.popen("gnuplot -persist","w") do |gp|
  gp.puts "set title 'General Boltzmann Machine'"
  gp.puts "set terminal postscript enhanced color eps"
  gp.puts "set xrange[0:1000]"
  gp.puts "set yrange[0:0.2]"
  gp.puts "set grid"
  gp.puts "set output 't-kl-s.eps'"
  gp.print "plot "
  gp.printf("'../RBM/dir-dat/t-kl.dat' title 'RBM P=2' w l,")
  gp.printf("'../GBM/dir-P2/t-kl.dat'  title 'GBM P=2' w l\n")
end

