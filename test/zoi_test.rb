require_relative 'dataset.rb'

class ZoiTest < Test::Unit::TestCase
  GENOME_PATH = "#{ROOT_PATH}/test/data/DNA_scaffolds_NODE_1.fa"
  attr_reader :dataset

  def setup
    @dataset = TestDataset.new
  end

  def create_zoi(sls: [])
    ContigElements::Zoi.new(DummyContig.new(sls), 100, 1000, '')
  end

  def test_sort_validation
    dataset.add_zoi([100, 150])
    zoi = dataset.contig.zoi.first
    zoi.validate
    zoi.check_defection

    assert_zoi_invalid zoi, :short_transcript
  end

  def test_no_hits_validation
    dataset.add_zoi([100, 500])
    zoi = dataset.contig.zoi.first
    zoi.validate
    zoi.check_defection

    assert_zoi_invalid zoi, :no_hits
  end

  def test_fused_genes_defection
    dataset.add_zoi([100, 500])
    dataset.add_hit([50, 150, 1], [400, 600, 1])

    zoi = dataset.contig.zoi.first
    zoi.validate
    zoi.check_defection

    assert_zoi_defective zoi, :fused_genes
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

  def test_hits_searching
    dataset.add_zoi([100, 200])
    dataset.add_hit([50, 110, 1], [120, 130, 1], [150, 250, 1])
    contig = dataset.contig
    zoi = contig.zoi.first

    assert_equal 2, zoi.blast_hits.count
    assert_equal [[120, 130], [150, 250]], zoi.blast_hits.map { |h| [h.start, h.finish] }
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
