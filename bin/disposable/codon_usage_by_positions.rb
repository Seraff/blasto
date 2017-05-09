#!/usr/bin/env ruby
require_relative '../lib.rb'

params = Slop.parse do |o|
  o.banner = "Get codon usage by positions, beginning from the right side"
  o.string '-fasta', '(required) fasta file'
  o.string '-gff', '(required) gff file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :gff, :fasta

puts 'Reading fasta file'
fasta = FastaReader.new(params[:fasta]).to_h
gff = Bio::GFF::GFF3.new(File.open(params[:gff]))

pb = ProgressBar.create title: 'Generating sequences',
                        starting_at: 0,
                        total: gff.records.count,
                        format: "%a %e %P% Processed: %c from %C"
usage = {}
extra_usage = {}

gff.records.each do |r|
  seq = fasta[r.seqname]
  raise "Node not found: #{r.seqname}" unless seq

  interval = seq[r.start-1..r.end-1].to_s.upcase
  if r.reverse_strand?
    interval = Bio::Sequence::NA.new(interval).reverse_complement.to_s.upcase
  end

  interval = interval.split('').each_slice(3).to_a.reverse.map(&:join)

  interval.each_with_index do |codon, i|
    usage[i] ||= {}

    usage[i][codon] ||= 0
    usage[i][codon] += 1

    usage[i][:total] ||= 0
    usage[i][:total] += 1
  end

  pb.increment
end

outfile = File.open "results/codon_usage_by_positions.json", 'w'
outfile.write usage.to_json
