#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Append "
  o.string '--fasta', 'Peptides fasta file, source for seqid'
  o.string '--gff', 'Annotation file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

peptides = {}
Bio::FlatFile.open(Bio::FastaFormat, params[:fasta]).each do |e|
  peptides[e.entry_id.to_s] = e.seq
end

outfile_path = append_to_filename(params[:gff], 'with_peptides')
outfile = File.open outfile_path, 'w'

File.open(params[:gff], 'r').each_line do |entry|
  if entry.start_with? '#'
    outfile.puts entry
  else
    comment = entry.split("\t").last
    seqid = comment.match(/ID=(?<seqid>\d+)_/)[:seqid]
    outfile.puts "#{entry.chomp}_#{peptides[seqid]}"
  end
end

outfile.close

