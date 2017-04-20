#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Convert blast hits file with amino acid coordinates to gff annotation with nucleotide indicies\n"
  o.string '-in', '(required) Input file.'
  o.string '-out', 'Output GFF file. If not provided, puts the output gff file to the input file folder.'
  o.string '-t', '(required) Source for generating gff file (query|subject)'
  o.bool '-m', '--merge', 'Merge close blast hits with the same query, subject and frame'
  o.integer '--max_distance', 'Max distance between close hits for merging'
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
outfile_path = params[:out] || change_path(params[:in], new_ext: 'gff')

reader = BlastReader.new input_file_path

if params[:m]
  merged_path = change_path(params[:in], append: "merged_#{SecureRandom.hex}")
  merged_file = File.open(merged_path, 'w')
  reader.merge_hits merged_file, target: params[:t], max_distance: params[:max_distance] || 256
  reader.close

  reader = BlastReader.new merged_path
end

gff_file = File.open outfile_path, 'w'
gff_file.puts "##gff-version 3"

pb = ProgressBar.create(title: 'Converting', starting_at: 0, total: reader.hits_count)

reader.each_hit do |hit|
  gff = hit.to_gff params[:t]
  gff_array = gff.split("\t")

  gff_file.puts gff_array.join("\t")
  pb.increment
end

puts "Finished"

gff_file.close
reader.file.close
File.delete(merged_file) if merged_file
