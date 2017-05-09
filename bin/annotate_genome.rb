#!/usr/bin/env ruby
require_relative 'lib.rb'

params = Slop.parse do |o|
  o.banner = "P57 genome annotator"
  o.string '-g', '--genome', '(required) Genome .fasta file.'
  o.string '-b', '--hits', '(required) Blast hits .csv file.'
  o.string '-r', '--reads', '(required) Reads .bam file.'
  o.bool '--skip-preparation', 'Use already prepared blast hits (if provided)'
  o.string '-o', 'Output gff file.'
  o.string '-t', '(required) Target for generating frame extracting (query|subject).'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :genome

unless params[:'skip-preparation']
  assure_params_provided params, :hits, :reads, :t
end

annotator = Annotator.new genome_path: params[:genome],
                          hits_path: params[:hits],
                          reads_path: params[:reads],
                          skip_preparation: params[:'skip-preparation']

annotator.prepare target: params[:t], mode: :genome # TODO: dynamically detect mode

# annotator.annotate

