#!/usr/bin/env ruby
ENV['BLASTO_ENV'] = 'TEST'

base_dir = File.expand_path(File.join(File.dirname(__FILE__), ".."))
test_dir = File.join(base_dir, "test")
bin_dir = File.join(base_dir, "bin")

require File.join(bin_dir, 'lib.rb')
require 'bundler/setup'
Bundler.require

exit Test::Unit::AutoRunner.run(true, test_dir)

