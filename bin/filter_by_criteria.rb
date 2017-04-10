#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Filters Blast Report report .csv file by evalue and qlen == qend"
  o.string '-in', '--input_file', '(required) Input csv file'
  o.string '-out', '--output_file', 'Output csv file'
  o.string '-t', '--threshold', '(required) Every hit with evalue >= threshold will be dropped'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

if !params[:in] || !params[:t]
  puts "Not all required options provided, use --help for more info"
  exit
end

def fits?(h)
  return true
  fits = h.data[:evalue].to_f < THRESHOLD
  # fits &&= h.data[:qlen].to_i == h.data[:qend].to_i
  # fits &&= h.data[:qstart].to_i == 1
  fits
end

THRESHOLD = params[:t].to_f

input_path = params[:in]
output_path = params[:out] || append_to_filename(input_path, 'filtered_by_evalue')

r = BlastReader.new input_path
pb = ProgressBar.create(title: 'Filtering by evalue', starting_at: 0, total: r.hits_count)

r.modify_hits(output_path) do |h|
  pb.increment
  next unless fits? h
  true
end
