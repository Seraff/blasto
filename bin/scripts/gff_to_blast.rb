#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Convert blast hits file with amino acid coordinates to gff annotation with nucleotide indicies\n"
  o.string '-in', '(required) Input file.'
  o.string '-out', 'Output GFF file.'
  o.string '-t', '--target', '(required) Source for generating gff file (query|subject)'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :in, :t

unless %w(query subject).include? params[:t]
  puts "Wrong source option, use --help for more info"
  exit 0
end

input_file_path = params[:in]
output_file_path = params[:out] || change_path(params[:in], new_ext: 'csv')

headers = {
           qseqid: nil,
           qlen: 0,
           sseqid: nil,
           slen: 0,
           length: 0,
           evalue: 0,
           pident: 0,
           bitscore: 0,
           mismatch: 0,
           gaps: 0,
           qstart: 0,
           qend: 0,
           sstart: 0,
           send: 0
          }

File.open(output_file_path, 'w') do |outfile|
  outfile.puts headers.keys.join(',')

  File.open(input_file_path, 'r').each do |line|
    next if line.start_with? '#'
    headers[:qseqid] = SecureRandom.hex
    headers[:sseqid] = SecureRandom.hex

    splitted = line.split("\t")

    seqid = splitted[0]
    start = splitted[3]
    finish = splitted[4]

    if params[:t] == 'query'
      headers[:qseqid] = seqid
      headers[:qstart] = start
      headers[:qend] = finish
    else
      headers[:sseqid] = seqid
      headers[:sstart] = start
      headers[:send] = finish
    end

    outfile.puts headers.values.join(',')
  end
end

puts 'Finished'
