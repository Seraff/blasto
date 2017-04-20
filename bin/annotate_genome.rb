#!/usr/bin/ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "P57 genome annotator"
  o.string '-g', '(required) Genome .fasta file.'
  o.string '-b', '(required) Blast hits .csv file.'
  o.string '-r', '(required) Reads .bam file.'
  o.string '-o', 'Output gff file.'
  o.string '-t', '(required) Source for generating frame extracting (query|subject).'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :g, :b, :r, :t

annotator = Annotator.new genome_path: params[:g],
                          hits_path: params[:b],
                          reads_path: params[:r]

annotator.prepare target: params[:t], mode: :genome # TODO: dynamically mode

annotator.annotate

