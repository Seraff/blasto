#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Filters Blast Report report .csv file by evalue"
  o.string '-gff', '(required) gff file with blast hits'
  o.string '-fasta', '(required) fasta file with genome'
  o.string '-out', 'Output blast hits gff file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

if !params[:gff] || !params[:fasta]
  puts "Not all required options provided, use --help for more info"
  exit
end

outpath = params[:out] || append_to_filename(params[:gff], 'before_stops')

fasta = {}
FastaReader.new(params[:fasta]).each do |contig|
  fasta[contig.entry_id] = contig.seq
end

out_file = File.open(outpath, 'w')
out_file.puts "##gff-version 3"

gff = Bio::GFF::GFF3.new(File.open(params[:gff]))
pb = ProgressBar.create(title: 'Choosing hits', starting_at: 0, total: gff.records.count)

gff.records.each do |r|
  pb.increment

  raise "Wrong fasta file: unfound sequence #{r.seqname}" unless fasta[r.seqname]
  seq = fasta[r.seqname]

  if [1, 2, 3].include? r.frame
    codon = seq[(r.end)..(r.end-1+3)]
  elsif [4, 5, 6].include? r.frame
    codon = seq[(r.start-1-3)..(r.start-1-1)]
    codon = Bio::Sequence::NA.new(codon).reverse_complement.upcase
  end

  next if codon.size != 3

  if %w(TAA TGA TAG).include? codon.upcase
    out_file.puts r.to_s
  end
end

out_file.close
