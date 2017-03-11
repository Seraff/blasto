#!/usr/bin/ruby
require_relative 'lib.rb'

MAX_DIST = 64

`bedtools cluster -d #{MAX_DIST} -i data/test/p57_DNA_hits_nucleotides_sorted.gff > data/test/p57_DNA_hits_nucleotides_sorted_clustered.gff`

outfile = File.open('data/test/p57_DNA_hits_nucleotides_sorted_clustered_merged.gff', 'w')

current_data = nil

File.open('data/test/p57_DNA_hits_nucleotides_sorted_clustered.gff', 'r').each do |entry|
  data = entry.split("\t")

  unless current_data
    current_data = data
    next
  end

  if current_data[-1] == data[-1]
    current_data[3] = data[3].to_i if data[3].to_i < current_data[3].to_i
    current_data[4] = data[4].to_i if data[4].to_i > current_data[4].to_i
  else
    outfile.puts current_data[0..-2].join("\t")
    current_data = data
  end
end

outfile.close
puts 'ok'
