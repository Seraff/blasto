#!/usr/bin/env ruby
require_relative '../lib.rb'

data = {}

File.open('results/annotator/annotation.gff').each do |l|
  s = l.split("\t")
  data[s[0]] ||= []
  data[s[0]] << [s[3].to_i, s[4].to_i]
end

problems = []

data.each do |contig, zois|
  zois = zois.sort_by { |e| e[0] }

  zois.each_with_index do |z, i|
    next if i == 0

    if [(z[0]..z[1]), (zois[i-1][0]..zois[i-1][1])].intersection
      problems << "#{contig} #{[(z[0]..z[1]), (zois[i-1][0]..zois[i-1][1])]}"
    end
  end
end

puts 'done'

puts problems.map { |p| "#{p}\n" }
