#!/usr/bin/env ruby
require_relative '../lib.rb'
require 'matplotlib/pyplot'

params = Slop.parse do |o|
  o.banner = "Visualize JSON file"
  o.string '-in', '(required) json file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

STOPS = ['TGA', 'TAA', 'TAG']
NON_STOPS = ['TGG', 'GAA']

def get_stats_by_usage(usage, indexes, codon)
  vals = [nil]

  indexes.each do |i|
    unless usage[i]
      vals << 0
      next
    end

    vals << (usage[i][codon].to_f/usage[i]['total'].to_f)*100
  end

  vals
end

f = File.open(params[:in], 'r')
usage = JSON.parse(f.read)
usage = usage.map { |k, v| [k.to_i, v] }.to_h

MAX_X = 100

plt = Matplotlib::Pyplot
plt.ylim(0, 10)
plots = []

(STOPS).each do |codon|
  vals = get_stats_by_usage usage, [*1..MAX_X], codon
  plots << plt.plot([*-MAX_X..0], vals.reverse, label: codon)
end

plt.legend
plt.show
