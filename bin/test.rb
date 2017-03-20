#!/usr/bin/ruby
require_relative 'lib.rb'

# generating fasta file from peptides cluster or something

hits_path = 'data/peptides/all_pepti_p57_genome_bl_report_best.csv'
identical_path = 'tmp/all_pepti_p57_genome_bl_report_identical.csv'

def hit_identical?(h)
  h.data[:pident].to_f == 100 &&
    h.data[:length].to_i == h.data[:qend].to_i &&
    h.data[:qstart].to_i == 1
end

r = BlastReader.new hits_path

r.modify_hits(identical_path) do |h|
  next unless hit_identical? h
  true
end

r = BlastReader.new hits_path

r.modify_hits('tmp/all_pepti_p57_genome_bl_report_not_identical.csv') do |h|
  next if hit_identical? h
  true
end

# convert to gff
gff_path = append_to_filename identical_path, 'nucleotides'
gff_path = [gff_path.split('.')[0..-2] + ['gff']].flatten.join('.')
`bin/translated_hits_to_gff.rb -m genome -t subject -in #{identical_path} -out #{gff_path}`

# sort
sorted_path = append_to_filename(gff_path, 'sorted')
`sort -n -k4,5 #{gff_path} > #{sorted_path}`

# cluster
clustered_path = append_to_filename(sorted_path, 'clustered_merged')
`bin/cluster_and_merge.rb -in #{sorted_path} -out #{clustered_path}`

# write results
clusters = {}
Bio::GFF.new(File.open(clustered_path)).records.each do |rec|
  next if rec.comment.to_s.start_with? '#'

  clusters[rec.seqname] ||= []
  clusters[rec.seqname] << [rec.start.to_i, rec.end.to_i]
end

outfile = File.open('tmp/results.fasta', 'w')

Bio::FlatFile.open(Bio::FastaFormat, 'data/DNA_scaffolds.fa').each do |contig|
  clusters[contig.entry_id].to_a.each do |cluster|
    outfile.puts "> #{contig.entry_id}_#{cluster.first}_#{cluster.last}"
    outfile.puts contig.seq[(cluster.first-1)..(cluster.last-1)]
  end
end

outfile.close

