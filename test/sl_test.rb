class DummyContig
  attr_accessor :sl_mappings

  def initialize(sl_data)
    sls = sl_data.map { |e| ContigElements::Sl.new self, e[0], e[1], { coverage: e[2] } }
    @sl_mappings = ContigElementCollection.new(sls)
  end

  def seq
    (0..1000).map { |e| %w(G, T, A, C).sample }.join
  end
end

class SlTest < Test::Unit::TestCase
  GENOME_PATH = "#{ROOT_PATH}/test/data/DNA_scaffolds_NODE_1.fa"

  def setup
  end

  def teardown
  end

  def create_zoi(sls: [])
    ContigElements::Zoi.new(DummyContig.new(sls), 100, 1000, '')
  end

  def test_selecting_internal
    zoi = create_zoi(sls: [[110, 111, 1], [200, 201, 1]])
    inner = zoi.sl_mappings.detect { |e| e.start == 200 }
    assert(zoi.sl_mapping == inner)
  end

  def test_selecting_by_coverage
    zoi = create_zoi(sls: [[110, 111, 2], [200, 201, 1]])
    covered = zoi.sl_mappings.detect { |e| e.coverage == 2 }
    puts zoi.sl_mapping.start
    assert(zoi.sl_mapping == covered)
  end
end
