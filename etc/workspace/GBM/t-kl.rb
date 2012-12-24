#!/usr/bin/ruby

dirnames = ["dir-cmp"]

IO.popen("gnuplot -persist","w") do |gp|
  gp.puts "set title 'General Boltzmann Machine'"
  gp.puts "set terminal postscript enhanced color eps"
  gp.puts "set xrange[0:1000]"
  gp.puts "set yrange[0:0.2]"
  gp.puts "set grid"
  gp.puts "set output '#{dirnames}/t-kl-s.eps'"
  gp.puts "plot 'dir-P1/t-kl.dat' w l title 'P=1', 'dir-P2/t-kl.dat' w l title 'P=2', 'dir-P3/t-kl.dat' w l title 'P=3', 'dir-P4/t-kl.dat' w l title 'P=4'"
end

