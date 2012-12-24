#!/usr/bin/ruby


fname = ARGV[0]
open(fname, 'r') do |fr|
  fr.gets
  puts fr.read
end


