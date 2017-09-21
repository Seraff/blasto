require 'rubygems'
require 'bundler/setup'
require 'json'
require 'securerandom'
require 'bigdecimal'

Bundler.require

Config.load_and_set_settings('config/annotator.yml')

require 'pathname'

ROOT_PATH = File.dirname Pathname.new(File.absolute_path(__FILE__)).parent

Dir["#{ROOT_PATH}/lib/helpers/*.rb"].each {|file| require file }

# Annotator stuff

require_relative '../lib/annotator.rb'
require_relative '../lib/preparer.rb'
require_relative '../lib/bad_transcripts_logger.rb'
