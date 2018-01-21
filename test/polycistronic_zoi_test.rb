require_relative 'dataset.rb'
require_relative 'annotation_test_helper.rb'

class PolycistronicZoiTest < Test::Unit::TestCase
  include AnnotationTestHelper
  attr_reader :dataset, :contig

  def setup
    Settings.annotator.transcript_min_size = 6
    Settings.annotator.zoi_hit_searching_inner_threshold = 3
  end

  def teardown
    Settings.annotator.transcript_min_size = 100
    Settings.annotator.zoi_hit_searching_inner_threshold = 15
  end

  def test_simple_polycistronic_split
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------[-------------------------]------ # ZOI
      --------[----]------------------------- # hit, frame 1
      ---------------------[----]------------ # hit, frame 2
      ---------------[]---------------------- # SL, coverage 1
    }

    build_dataset_from_string(data)

    assert_true zoi.polycistronic?
    assert_equal [16], zoi.polycistronic_cutting_places
    assert_equal 2, zoi.blast_hit_groups.count
    assert_equal 2, zoi.blast_hits.count
  end

  def test_simple_polycistronic_split_without_sl
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------[-------------------------]------ # ZOI
      --------[----]------------------------- # hit, frame 1
      ---------------------[----]------------ # hit, frame 2
    }

    build_dataset_from_string(data)

    assert_true zoi.polycistronic?
    assert_equal [18], zoi.polycistronic_cutting_places
    assert_equal 2, zoi.blast_hit_groups.count
    assert_equal 2, zoi.blast_hits.count
  end

  def test_simple_polycistronic_split_same_frame_hits
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------[-------------------------]------ # ZOI
      --------[----]------------------------- # hit, frame 1
      --------------------[----]------------- # hit, frame 1
      ---------------[]---------------------- # SL, coverage 1
    }

    build_dataset_from_string(data)

    assert_true zoi.polycistronic?
    assert_equal [16], zoi.polycistronic_cutting_places
    assert_equal 2, zoi.blast_hit_groups.count
    assert_equal 2, zoi.blast_hits.count
  end

  def test_not_polycistronic_hits
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------[-------------------------]------ # ZOI
      --------[----]------------------------- # hit, frame 1
      --------------------[----]------------- # hit, frame 1
    }

    build_dataset_from_string(data)

    assert_false zoi.polycistronic?
    assert_equal [], zoi.polycistronic_cutting_places
    assert_equal 1, zoi.blast_hit_groups.count
    assert_equal 2, zoi.blast_hits.count

    assert_equal 9, zoi.blast_hit_begin
    assert_equal 26, zoi.blast_hit_end
  end
end
