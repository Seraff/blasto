require_relative 'dataset.rb'

class SlTest < Test::Unit::TestCase
  GENOME_PATH = "#{ROOT_PATH}/test/data/DNA_scaffolds_NODE_1.fa"
  attr_reader :dataset

  def setup
    @dataset = TestDataset.new
  end

  def teardown
  end

  def create_zoi(sls: [])
    ContigElements::Zoi.new(DummyContig.new(sls), 100, 1000, '')
  end

  def test_selecting_internal
    self.dataset.add_zoi([100, 1000])
    self.dataset.add_sl([110, 111, 1], [200, 201, 1])
    zoi = self.dataset.contig.zoi.first

    inner = zoi.sl_mappings.detect { |e| e.start == 200 }
    assert(zoi.sl_mapping == inner)
  end

  def test_selecting_by_coverage
    self.dataset.add_zoi([100, 1000])
    self.dataset.add_sl([110, 111, 1], [200, 201, 2])
    zoi = self.dataset.contig.zoi.first

    covered = zoi.sl_mappings.detect { |e| e.coverage == 2 }
    assert(zoi.sl_mapping == covered)
  end
end
