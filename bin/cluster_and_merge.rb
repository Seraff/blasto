#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Cluster and merge .gff file"
  o.string '-in', '--input_file', 'Input gff file'
  o.string '-out', '--output_file', 'Output gff file'
  o.string '-m', '--max_distance', 'Max distance between sequences for merging'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

unless params[:in]
  puts "-in option is required, use --help for more info"
  exit
end

MAX_DIST = params[:m] || 64

# sort
sorted_path = append_to_filename(params[:in], 'sorted')
`sort -n -k4,5 #{params[:in]} > #{sorted_path}`

clustered_path = append_to_filename(sorted_path, 'clustered')
`bedtools cluster -s -d #{MAX_DIST} -i #{sorted_path} > #{clustered_path}`

outpath = params[:out] || append_to_filename(params[:in], 'clusters')
outfile = File.open(outpath, 'w')
current_data = nil

File.open(clustered_path, 'r').each do |entry|
  data = entry.split("\t")

  unless current_data
    current_data = data
    next
  end

  if current_data[-1] == data[-1]
    current_data[3] = data[3].to_i if data[3].to_i < current_data[3].to_i
    current_data[4] = data[4].to_i if data[4].to_i > current_data[4].to_i
  else
    outfile.puts current_data[0..-2].join("\t")
    current_data = data
  end
end

`rm #{sorted_path} #{clustered_path}`
outfile.close
