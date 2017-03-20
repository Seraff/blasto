#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Append "
  o.string '-fasta', '(required) fasta file'
  o.string '-gff', '(required) gff file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :gff, :fasta

fasta = FastaReader.new(params[:fasta]).to_h
gff = Bio::GFF::GFF3.new(File.open(params[:gff]))
pb = ProgressBar.create(title: 'Gethering codon usage', starting_at: 0, total: gff.records.count)

total_usage = {}
real_stops_usage = {}

gff.records.each do |r|
  seq = fasta[r.seqname]
  raise "Node not found: #{r.seqname}" unless seq

  interval = seq[r.start-1..r.end-1]

  if r.reverse_strand?
    interval = Bio::Sequence::NA.new(interval).reverse_complement.upcase
  end

  Bio::Sequence::NA.new(interval).codon_usage.each do |c, count|
    aa = Bio::Sequence::NA.new(c).translate
    total_usage[aa] ||= {}
    total_usage[aa][c.upcase] ||= 0
    total_usage[aa][c.upcase] += count
  end

  # real stops
  if r.reverse_strand?
    stop = seq[(r.start-1-3)..(r.start-1-1)]
    stop = Bio::Sequence::NA.new(stop).reverse_complement.upcase
  else
    stop = seq[(r.end)..(r.end-1+3)]
  end

  real_stops_usage[stop] ||= 0
  real_stops_usage[stop] += 1

  pb.increment
end

puts "Total usage"
total_usage.each do |codon, hsh|
  puts codon
  hsh.each do |na, count|
    puts "\t#{na}: #{count}"
  end
end

puts "\nReal stops usage"
real_stops_usage.each do |na, count|
  puts "\t#{na}: #{count}"
end
