#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Tests annotated gener for connsistency (start/stop codons, length)"
  o.string '--genome', '(required) Genome fasta file'
  o.string '--gff', '(required) Annotation'
  o.on '-h', '--help', 'Print options' do
    puts o
  end
end

def check_codons(interval, strand)
  frame = strand == '+' ? 1 : 4
  translated = Bio::Sequence::NA.new(interval).translate(frame).to_s

  return false if translated[0] != 'M' || translated[-1] != '*'
  true
end

def check_multiplicity(interval)
  return false if interval.size % 3 != 0
  return true
end

def check_length(interval)
  return false if interval.size < Settings.annotator.gene_min_size
  true
end

fasta = FastaReader.new(params[:genome]).to_h
gff = Bio::GFF::GFF3.new(File.open(params[:gff]))
pb = ProgressBar.create title: 'Checking sequences',
                        starting_at: 0,
                        total: gff.records.count

errors_cnt = 0

gff.records.each do |r|
  seq = fasta[r.seqname]
  start = r.start-1
  finish = r.end-1
  start -= 3 if r.strand == '-'
  finish += 3 if r.strand == '+'
  interval = seq[start..finish].to_s.upcase

  unless check_codons interval, r.strand
    puts "#{r.seqname} (#{r.start}-#{r.end}): interval not valid (#{interval})"
    errors_cnt += 1
  end

  unless check_multiplicity interval
    puts "#{r.seqname} (#{r.start}-#{r.end}): interval is not a multiple of three"
    errors_cnt += 1
  end

  unless check_length interval
    puts "#{r.seqname} (#{r.start}-#{r.end}): interval is short"
    errors_cnt += 1
  end

  pb.increment
end

puts "#{errors_cnt} errors found"
