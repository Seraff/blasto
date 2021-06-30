require_relative 'dataset.rb'
require_relative 'annotation_test_helper.rb'

class AnnotationTest < Test::Unit::TestCase
  include AnnotationTestHelper
  attr_reader :dataset, :contig

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
    zoi.validate
    zoi.check_defection

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
    zoi.validate
    zoi.check_defection

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

    result = %q{
      ---------[-------------]--------------- # annotated gene
    }

    build_dataset_from_string data
    zoi.validate
    zoi.check_defection

    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.annotate
    assert_false zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.validation_errors.include?(:cannot_detect_start)
    assert_annotation result, direction: '+'
  end

  def test_stop_detection_fail
    data = %q{
      TAGGTACATATGGGTACATGCATCTAGGATGTATAAGGT # genome
      ---------[M]---------------------[*]---
      -------[-------------------]----------- # ZOI
      -------------[----------]-------------- # hit, frame 1
      ----[]--------------------------------- # SL, coverage 1
    }

    result = %q{
      ---------[----------------]------------ # annotated gene
    }

    build_dataset_from_string data
    zoi.validate
    zoi.check_defection

    assert_true zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.annotate
    assert_false zoi.valid?
    assert_false zoi.defective?
    assert_true zoi.validation_errors.include?(:cannot_detect_stop)
    assert_annotation result, direction: '+'
  end

  def test_gene_begin_detection_fail_reverse
    data = %q{
      TAGGTACTTAAGGGTACATGCATCTGCTTGGATGTAGGT # genome
      -------[*]-----------------------------
      -----[----------------------------]---- # ZOI
      -------------[----------------]-------- # hit, frame 6
    }

    result = %q{
      ----------[----------------------]----- # annotated gene
    }

    build_dataset_from_string data
    zoi.validate
    zoi.check_defection

    assert_true zoi.valid?
    assert_true zoi.annotate
    assert_false zoi.valid?
    assert_true zoi.validation_errors.include?(:cannot_detect_start)
    assert_annotation result, direction: '-'
  end

  def test_gene_end_detection_fail_reverse
    data = %q{
      TAGGTACATCAGGGTACATGCATCTGGGCATATGTAGGT # genome
      ----------------------------[M]--------
      -----[----------------------------]---- # ZOI
      -------------[----------]-------------- # hit, frame 6
    }

    result = %q{
      -------[----------------------]-------- # annotated gene
    }

    build_dataset_from_string data
    zoi.validate
    zoi.check_defection

    assert_true zoi.valid?
    assert_true zoi.annotate
    assert_false zoi.valid?
    assert_true zoi.validation_errors.include?(:cannot_detect_stop)
    assert_annotation result, direction: '-'
  end

  def test_correct_sl_absence
    data = %q{
      TAGGTACATATGGGTACATGCATCTAATAGGATGTAGGTGTAGGT # genome
      ---------[M]------------[*]------------------
      -------[-------------------]----------------- # ZOI
      -------------[----------]-------------------- # hit, frame 1
      -------------------------------------------[] # SL, coverage 3
    }

    build_dataset_from_string data
    zoi.validate
    zoi.check_defection

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
    zoi.validate
    zoi.check_defection

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
    zoi.validate
    zoi.check_defection

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
    zoi.validate
    zoi.check_defection

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
    zoi.validate
    zoi.check_defection

    assert_false zoi.begin_torn?
    assert_true zoi.end_torn?

    assert_true zoi.annotate

    assert_true zoi.valid?
    assert_true zoi.defective?
    assert_annotation result, direction: '-'
  end
end
