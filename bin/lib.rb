require 'rubygems'
require 'bundler/setup'
require 'json'
require 'securerandom'

Bundler.require

Config.load_and_set_settings('config/annotator.yml')

require 'pathname'

ROOT_PATH = File.dirname Pathname.new(File.absolute_path(__FILE__)).parent

require_relative '../lib/helpers/helpers.rb'
require_relative '../lib/helpers/blast_reader.rb'
require_relative '../lib/helpers/fasta_reader.rb'
require_relative '../lib/helpers/reads_statistics.rb'
require_relative '../lib/helpers/gff_clusterizer.rb'

# Annotator stuff

require_relative '../lib/annotator.rb'
require_relative '../lib/preparer.rb'
