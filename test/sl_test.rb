require_relative 'dataset.rb'

class SlTest < Test::Unit::TestCase
  GENOME_PATH = "#{ROOT_PATH}/test/data/DNA_scaffolds_NODE_1.fa"
  attr_reader :dataset

  def setup
    @dataset = TestDataset.new
  end

  def teardown
  end

  #TODO
end
