require_relative 'dataset.rb'

class AnnotationTest < Test::Unit::TestCase
  include AnnotationTestHelper
  attr_reader :dataset, :contig

  def setup
    Settings.annotator.transcriptome_min_size = 6
  end

  def teardown
    Settings.annotator.transcriptome_min_size = 100
  end

  def test_forward_frame_correct_annotation
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ---------[M]------------[*]------------
      -------[-------------------]----------- # ZOI
      -------------[----------]-------------- # hit, frame 1
      ----[]--------------------------------- # SL, coverage 1
    }

    result = %q{
      ---------[-------------]--------------- # annotated gene
    }

    build_dataset_from_string(data)
    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.annotate
    assert_annotation result, direction: '+'
  end

  def test_reverse_frame_correct_annotation
    data = %q{
      TAGGTACATTAGGGTACATGCATCTGCATGGATGTAGGT # genome
      --------[*]---------------[M]----------
      -[------------------------------------] # ZOI
      --------------[----------]------------- # hit, frame 5
      ------------------------------------[]- # SL, coverage 1
    }

    result = %q{
      -----------[----------------]---------- # result
    }

    build_dataset_from_string(data)
    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.annotate
    assert_annotation result, direction: '-'
  end

  def test_start_detection_fail
    data = %q{
      TAGGTACATAAGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------------------------[*]------------
      -------[-------------------]----------- # ZOI
      -------------[----------]-------------- # hit, frame 1
      ----[]--------------------------------- # SL, coverage 1
    }

    build_dataset_from_string data

    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_false zoi.annotate
    assert_false zoi.valid?
    assert_equal :cannot_detect_start, zoi.validation_error
    assert_false zoi.defective?
  end

  def test_stop_detection_fail
    data = %q{
      TAGGTACATATGGGTACATGCATCTAGGATGTATAAGGT # genome
      ---------[M]---------------------[*]---
      -------[-------------------]----------- # ZOI
      -------------[----------]-------------- # hit, frame 1
      ----[]--------------------------------- # SL, coverage 1
    }

    build_dataset_from_string data

    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_false zoi.annotate
    assert_false zoi.valid?
    assert_equal :cannot_detect_stop, zoi.validation_error
    assert_false zoi.defective?
  end

  def test_correct_sl_absence
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGT # genome
      ---------[M]------------[*]------------
      -------[-------------------]----------- # ZOI
      -------------[----------]-------------- # hit, frame 1
      -------------------------------------[] # SL, coverage 1
    }

    build_dataset_from_string data

    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.annotate
    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_equal :blast_hit_has_another_frame, zoi.defection_reason
  end

  protected

  def build_dataset(seq)
    @dataset = TestDataset.new(seq)
    @contig = @dataset.contig
  end

  def zoi
    @zoi ||= dataset.contig.zoi.first
  end
end
