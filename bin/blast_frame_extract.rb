#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "Extract frame from qseqid/sseqid to the sparate column"
  o.string '-in', '(required) Input file.'
  o.string '-out', 'Output file. If not provided, puts the output blast hits file to the input file folder.'
  o.string '-t', '(required) Source for generating frame extracting (query|subject).'
  o.string '-d', 'Delimiter of blast hits file'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :in, :t

delimiter = params[:d] || ','

out_path = params[:out] || change_path(params[:in], append: 'frames_extracted')
out_file = File.open(out_path, 'w')

reader = BlastReader.new(params[:in])

new_column_name = case params[:t]
when 'query'
  'qframe'
when 'subject'
  'sframe'
end

out_file.puts "#{reader.headers.join(delimiter)}#{delimiter}#{new_column_name}"

pb = ProgressBar.create(title: 'Extracting frame', starting_at: 0, total: reader.hits_count)

reader.each_hit do |hit|
  frame = hit.detect_frame(params[:t])

  seqid_key = BlastHit::TARGET_KEYS[params[:t].to_sym][:id]
  hit.data[seqid_key].gsub!(/\_\d+\z/, '')

  out_file.puts [hit.to_csv(delimiter: delimiter), frame].join(delimiter)
  pb.increment
end
