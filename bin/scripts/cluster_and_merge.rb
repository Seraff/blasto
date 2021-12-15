#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Cluster and merge .gff file"
  o.string '-in', '--input_file', '(required) Input gff file'
  o.string '-out', '--output_file', 'Output gff file'
  o.string '-m', '--max_distance', 'Max distance between sequences for merging'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :in
assure_file_param_has_extension params, :in, :gff
assure_file_param_has_extension params, :out, :gff if params[:out]

max_dist = params[:m] || 64
outpath = params[:out] || append_to_filename(params[:in], 'clusters')

clusterizer = Clusterizers::Gff.new input: params[:in],
                                    output: outpath,
                                    max_distance: max_dist

clusterizer.perform

puts 'Finished'
