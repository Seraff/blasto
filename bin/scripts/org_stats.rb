#!/usr/bin/env ruby
require_relative '../lib.rb'

params = Slop.parse do |o|
  o.banner = "Get per cluster organism stats in prepared annotation"
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

TARGET = :query

def orgs_of_cluster(cluster)
  cluster.data.map { |h| h.data.organism(TARGET) }.compact
end

annotator = Annotator.new

cnt = 0
bad_cnt = 0

pb = ProgressBar.create title: 'Calculating',
                        starting_at: 0,
                        total: annotator.count,
                        format: "%a %e %P% Processed: %c from %C"

annotator.each_contig do |cont|
  cont.blast_hit_clusters.each do |clust|
    orgs = orgs_of_cluster clust
    if orgs.group_by { |e| e }.values.map(&:count).max.to_i > 1
      cnt += 1
    else
      bad_cnt += 1
    end
  end
  pb.increment
end

puts "Good: #{cnt}"
puts "Bad: #{bad_cnt}"

