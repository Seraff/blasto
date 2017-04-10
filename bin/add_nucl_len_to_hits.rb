#!/usr/bin/ruby
require_relative 'lib.rb'

nucl_fasta_path = select_file(subject: "nucleotide FASTA file")
bh_path = select_file(subject: "Blast hits without nucleotides length")
output_path = append_to_filename bh_path, 'with_nlen'

fasta_data = Bio::FastaFormat.open(nucl_fasta_path)
                             .map{|e| [e.first_name, e.seq.size] }.to_h

outfile = File.open output_path, 'w'
bh_file = File.open(bh_path, 'r')

bh_file.each do |line|
  data = line.split(',')
  if data.first == 'qseqid'
    outfile.puts data.join(',')
    next
  end

  qseqid = data.first(3).last
  formatted_qseqid = qseqid.split('_')[0..-2].join('_')
  nlen = fasta_data[formatted_qseqid]
  raise "Cannot detect nlen for #{formatted_qseqid}" unless nlen

  data[2] = [qseqid.split('_')[0..-2], "length_#{nlen}", qseqid.split('_').last].flatten.join('_')

  outfile.puts data.join(',')
end

puts "Finished"

outfile.close
bh_file.close
