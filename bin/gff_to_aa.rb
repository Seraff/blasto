#!/usr/bin/env ruby
require 'slop'
require 'bio'
require 'pry'
require 'securerandom'
require 'ruby-progressbar'

params = Slop.parse do |o|
  o.banner = "Translate .gff annotation to amino acid .fasta"
  o.string '-a', '--annotation', '(required) gff annotation file'
  o.string '-g', '--genome', '(required) genome fasta file'
  o.string '-o', '--output', 'output fasta file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

outfile_path = params[:output] || params[:annotation].gsub('.gff', '_peptides.fasta')

genome = {}

Bio::FlatFile.new(Bio::FastaFormat, File.open(params[:genome])).each do |e|
  genome[e.entry_id] = e
end

annotation = Bio::GFF.new File.open(params[:annotation])
pb = ProgressBar.create(title: 'Translating', starting_at: 0, total: annotation.records.count)

table = Bio::CodonTable[1]
table['taa'] = 'E'
table['tag'] = 'E'
table['tga'] = 'W'

File.open(outfile_path, 'w') do |out_f|
  annotation.records.each do |gff|
    entry = Bio::Sequence.auto(genome[gff.seqname].seq)[(gff.start.to_i-1)..(gff.end.to_i-1)]
    frame = gff.strand == '+' ? 1 : 4

    entry = entry.translate(frame)

    out_f.puts ">#{SecureRandom.hex}|#{gff.seqname}|#{gff.start}-#{gff.end}|#{gff.attributes.keys.join(',')}"
    out_f.puts entry

    pb.increment
  end
end
