require_relative 'dataset.rb'
require_relative 'annotation_test_helper.rb'

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
    assert_true zoi.validation_errors.include?(:cannot_detect_start)
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
    assert_true zoi.validation_errors.include?(:cannot_detect_stop)
    assert_false zoi.defective?
  end

  def test_correct_sl_absence
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGTGTAGGT # genome
      ---------[M]------------[*]------------------
      -------[-------------------]----------------- # ZOI
      -------------[----------]-------------------- # hit, frame 1
      ---------------------------------------------[] # SL, coverage 1
    }

    build_dataset_from_string data

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_true zoi.annotate
    assert_true zoi.valid?
    assert_true zoi.defective?

    assert_true zoi.defection_reasons.include?(:hit_in_another_frame)
  end

  ## Torn cases

  def test_begin_torn_forward_zoi
    data = %q{
      TAGGTACATAGGGGTACATGCATCTAATAGGATGTAGGT # genome
      ------------------------[*]------------
      -------[-------------------]----------- # ZOI
      [----------------------]--------------- # hit, frame 1
      ----[]--------------------------------- # SL, coverage 1
    }

    result = %q{
      [----------------------]--------------- # annotated gene
    }

    build_dataset_from_string data

    assert_true zoi.begin_torn?
    assert_false zoi.end_torn?

    assert_true zoi.annotate

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_annotation result, direction: '+'
  end

  def test_begin_torn_reverse_zoi
    data = %q{
      TAGGTACAGTTAGGGTACATGCATCTGCCTGGATGTAGG # genome
      ---------[*]---------------------------
      -------[-------------------]----------- # ZOI
      ---------------[----------------------] # hit, frame 4
      ----[]--------------------------------- # SL, coverage 4
    }

    result = %q{
      ------------[-------------------------] # annotated gene
    }

    build_dataset_from_string data

    assert_true zoi.begin_torn?
    assert_false zoi.end_torn?

    assert_true zoi.annotate

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_annotation result, direction: '-'
  end

  def test_end_torn_forward_zoi
    data = %q{
      TAGGTACATATGGGTACATGCATCTСATAGGATGTAGGT # genome
      ---------[M]---------------------------
      -------[-------------------]----------- # ZOI
      -----[--------------------------------] # hit, frame 1
      --------------------------------------- # SL, coverage 1
    }

    result = %q{
      ---------[----------------------------] # annotated gene
    }

    build_dataset_from_string data

    assert_false zoi.begin_torn?
    assert_true zoi.end_torn?

    assert_true zoi.annotate

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_annotation result, direction: '+'
  end

  def test_end_torn_reverse_zoi
    data = %q{
      TAGGTACATGGTACATGCATCCATTСATAGGATGTAGGT # genome
      -----------------[M]--------------------
      -------[------------------]------------- # ZOI
      [-------------]------------------------- # hit, frame 5
      ---------------------------------------- # SL, coverage 5
    }

    result = %q{
      [------------------]-------------------- # annotated gene
    }

    build_dataset_from_string data

    assert_false zoi.begin_torn?
    assert_true zoi.end_torn?

    assert_true zoi.annotate

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_annotation result, direction: '-'
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
