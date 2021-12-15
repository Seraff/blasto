#!/usr/bin/env ruby
require_relative '../lib.rb'

annotator = Annotator.new

pb = ProgressBar.create title: 'Calculating',
                        starting_at: 0,
                        total: annotator.count,
                        format: "%a %e %P% Processed: %c from %C"

total_count = 0
counter = 0

annotator.each_contig do |contig|
  contig.zoi.each do |z|
    other_hits = z.blast_hits.select do |bh|
      bh.data.data[:qseqid] != z.best_blast_hit.data.data[:qseqid]
    end

    best_org = z.best_blast_hit.organism
    other_orgs = other_hits.map(&:organism)

    raise 'error' if other_orgs.count == z.blast_hits.count

    counter += 1 if other_orgs.include?(best_org)
    total_count += 1
  end

  pb.increment
end

print "Found: #{counter}/#{total_count}"
