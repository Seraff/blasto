#!/usr/bin/ruby
require_relative '../lib.rb'

params = Slop.parse do |o|
  o.string '-genome', '(required) Genome fasta file'
  o.string '-hits', '(required) Input csv file'
  o.string '-r', 'Recalculate'
  o.on '-h', '--help', 'Print options' do
    puts o
    exit
  end
end

assure_params_provided params, :genome, :hits
recalc = params[:r] ? true : false

# filter by evalue and qlen&qend
filtered_path = change_path(params[:hits], new_dir: 'tmp', append: 'filtered')

if !File.exists?(filtered_path) || recalc
  system "bin/filter_by_criteria.rb -t 0.0001 -in #{params[:hits]} -out #{filtered_path}"
end

# back translate to gff
non_tr_path = change_path(params[:hits], new_dir: 'tmp', append: 'non_translated', new_ext: 'gff')

if !File.exists?(non_tr_path) || recalc
  system "bin/translated_hits_to_gff.rb -m genome -t subject -in #{filtered_path} -out #{non_tr_path}"
end

# # choose before stops
# bs_path = change_path(params[:hits], new_dir: 'tmp', append: 'before_stops', new_ext: 'gff')

# if !File.exists?(bs_path) || recalc
#   system "bin/choose_bh_before_stops.rb -fasta #{params[:genome]} -gff #{non_tr_path} -out #{bs_path}"
# end

# cluster and merge
clusters_path = change_path(params[:hits], new_dir: 'tmp', append: 'clusters', new_ext: 'gff')

if !File.exists?(clusters_path) || recalc
  system "bin/cluster_and_merge.rb -m 0 -in #{non_tr_path} -out #{clusters_path}"
end

system "bin/codon_usage.rb -fasta #{params[:genome]} -gff #{clusters_path}"
