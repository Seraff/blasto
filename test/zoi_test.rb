require_relative 'dataset.rb'

class ZoiTest < Test::Unit::TestCase
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

  def test_sort_validation
    dataset.add_zoi([100, 150])
    zoi = dataset.contig.zoi.first

    assert_zoi_invalid zoi, :short
  end

  def test_no_hits_validation
    dataset.add_zoi([100, 500])
    zoi = dataset.contig.zoi.first

    assert_zoi_invalid zoi, :no_hits
  end

  def test_fused_genes_defection
    dataset.add_zoi([100, 500])
    dataset.add_hit([50, 120, 1], [450, 600, 1])

    zoi = dataset.contig.zoi.first

    assert_zoi_defective zoi, :fused_genes
  end

  def test_hits_inconsistency_exception
    dataset.add_zoi([100, 500])
    dataset.add_hit([50, 120, 1], [450, 600, 2])

    assert_raise ContigElements::Zoi::ZoiHitsException do
      dataset.contig.zoi.first.valid?
    end
  end

  def test_sl_sorting
    dataset.add_zoi([100, 500])
    dataset.add_sl([90, 91, 1], [110, 111, 3], [480, 481, 2])
    contig = dataset.contig

    source_sls = contig.sl_mappings
    sorted_sls_left = contig.zoi.first.left_sls_sorted
    sorted_sls_right = contig.zoi.first.right_sls_sorted

    assert_equal sorted_sls_left[0], source_sls[1]
    assert_equal sorted_sls_left[1], source_sls[0]

    assert_equal sorted_sls_right[0], source_sls[2]
  end

  protected

  def assert_zoi_invalid(zoi, reason)
    assert zoi.invalid?
    assert_true zoi.validation_errors.include?(reason.to_sym)
  end

  def assert_zoi_defective(zoi, reason)
    assert zoi.defective?
    assert_true zoi.defection_reasons.include?(reason.to_sym)
  end
end
