#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "P57 genome annotator"
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

required_params = [:genome]
required_params += [:genome_reads, :blast_hits, :transcriptome, :blast_hit_target] unless Settings.annotator.skip_preparation

required_params.each do |name|
  raise "#{name} param is not provided in annotator.yml" unless Settings.annotator.send(name)
end

annotator = Annotator.new
annotator.prepare
annotator.annotate

