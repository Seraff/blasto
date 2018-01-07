require_relative 'dataset.rb'

class ContigSubsequenceTest < Test::Unit::TestCase
  attr_reader :dataset, :contig, :subseq

  def test_correct_subseq_seq
    build_dataset('GTAGTAGTA')
    build_subseq(1, 3)
    assert_equal subseq.seq.downcase, 'gta'
  end

  def test_correct_subseq_borders
    build_dataset('GTAGTAGTA')
    build_subseq(2, 5)
    assert_equal subseq.left_id, 2
    assert_equal subseq.right_id, 5
  end

  def test_custom_codon_table
    seq = 'TAG TGA TAA'
    build_dataset(seq.gsub(' ', ''))
    build_subseq(1, 9)

    aa = subseq.aa_seq_in_frame(1)
    assert_equal aa, '^^*'
  end

  def test_translation
    seq = 'GTACAGTCGGNNNTCATA'

    build_dataset(seq)
    build_subseq(1, seq.size)

    # GTA CAG TCG GNN NTC ATA
    assert_equal 'VQSXXI', subseq.aa_seq_in_frame(1)

    # G TAC AGT CGG NNN TCA TA
    assert_equal 'YSRXS', subseq.aa_seq_in_frame(2)

    # GT ACA GTC GGN NNT CAT A
    assert_equal 'TVXXH', subseq.aa_seq_in_frame(3)

    # TAT GAN NNC CGA CTG TAC
    assert_equal 'YXXRLY', subseq.aa_seq_in_frame(4)

    # T ATG ANN NCC GAC TGT AC
    assert_equal 'MXXDC', subseq.aa_seq_in_frame(5)

    # TA TGA NNN CCG ACT GTA C
    assert_equal '^XPTV', subseq.aa_seq_in_frame(6)
  end

  def test_subseq_frame_by_contig_frame
    seq = 'GTACAGTCGGNNNTCATAGACTAG'
    build_dataset(seq)

    ## Forward:
    # subsequence in 1 frame
    build_subseq(4, 6)
    assert_equal 1, @subseq.forward_subseq_frame
    assert_equal 1, @subseq.subseq_frame_by_contig_frame(1)
    assert_equal 2, @subseq.subseq_frame_by_contig_frame(2)
    assert_equal 3, @subseq.subseq_frame_by_contig_frame(3)

    # subsequence in 2 frame
    build_subseq(5, 7)
    assert_equal 2, @subseq.forward_subseq_frame
    assert_equal 3, @subseq.subseq_frame_by_contig_frame(1)
    assert_equal 1, @subseq.subseq_frame_by_contig_frame(2)
    assert_equal 2, @subseq.subseq_frame_by_contig_frame(3)

    # subsequence in 3 frame
    build_subseq(6, 8)
    assert_equal 3, @subseq.forward_subseq_frame
    assert_equal 2, @subseq.subseq_frame_by_contig_frame(1)
    assert_equal 3, @subseq.subseq_frame_by_contig_frame(2)
    assert_equal 1, @subseq.subseq_frame_by_contig_frame(3)

    ## Reverse:
    # subsequence in 4 frame
    build_subseq(22, 24)
    assert_equal 4, @subseq.reverse_subseq_frame
    assert_equal 4, @subseq.subseq_frame_by_contig_frame(4)
    assert_equal 5, @subseq.subseq_frame_by_contig_frame(5)
    assert_equal 6, @subseq.subseq_frame_by_contig_frame(6)

    # subsequence in 5 frame
    build_subseq(21, 23)
    assert_equal 5, @subseq.reverse_subseq_frame
    assert_equal 6, @subseq.subseq_frame_by_contig_frame(4)
    assert_equal 4, @subseq.subseq_frame_by_contig_frame(5)
    assert_equal 5, @subseq.subseq_frame_by_contig_frame(6)

    # subsequence in 6 frame
    build_subseq(20, 22)
    assert_equal 6, @subseq.reverse_subseq_frame
    assert_equal 5, @subseq.subseq_frame_by_contig_frame(4)
    assert_equal 6, @subseq.subseq_frame_by_contig_frame(5)
    assert_equal 4, @subseq.subseq_frame_by_contig_frame(6)
  end

  def test_aa_seq_in_contig_frame
    seq = 'GTACAGTCGGNNNTCATAGACTAG' # rev.c.: CTAGTCTATGANNNCCGACTGTAC
    build_dataset(seq)

    build_subseq(4, 12) # CAGTCGGNN (1/4 frame), rev.c: NNCCGACTG
    assert_equal 'QSX', @subseq.aa_seq_in_contig_frame(1) # GTA [CAG TCG GNN] NTC ATA GAC TAG
    assert_equal 'SR', @subseq.aa_seq_in_contig_frame(2) # G TA[C AGT CGG NN]N TCA TAG ACT AG
    assert_equal 'VX', @subseq.aa_seq_in_contig_frame(3) # GT A[CA GTC GGN N]NT CAT AGA CTA G
    assert_equal 'XRL', @subseq.aa_seq_in_contig_frame(4) # CTA GTC TAT GAN [NNC CGA CTG] TAC
    assert_equal 'XD', @subseq.aa_seq_in_contig_frame(5) # C TAG TCT ATG AN[N NCC GAC TG]T AC
    assert_equal 'PT', @subseq.aa_seq_in_contig_frame(6) # CT AGT CTA TGA N[NN CCG ACT G]TA C

    build_subseq(2, 10) # TACAGTCGG (2/5 frame), rev.c: CCGACTGTA
    assert_equal 'QS', @subseq.aa_seq_in_contig_frame(1) # G[TA CAG TCG G]NN NTC ATA GAC TAG
    assert_equal 'YSR', @subseq.aa_seq_in_contig_frame(2) # G [TAC AGT CGG] NNN TCA TAG ACT AG
    assert_equal 'TV', @subseq.aa_seq_in_contig_frame(3) # G[T ACA GTC GG]N NNT CAT AGA CTA G
    assert_equal 'RL', @subseq.aa_seq_in_contig_frame(4) # CTA GTC TAT GAN NN[C CGA CTG TA]C
    assert_equal 'DC', @subseq.aa_seq_in_contig_frame(5) # C TAG TCT ATG ANN N[CC GAC TGT A]C
    assert_equal 'PTV', @subseq.aa_seq_in_contig_frame(6) # CT AGT CTA TGA NNN [CCG ACT GTA] C

    build_subseq(3, 11) # ACAGTCGGN (3/6 frame), rev.c: NCCGACTGT
    assert_equal 'QS', @subseq.aa_seq_in_contig_frame(1) # GT[A CAG TCG GN]N NTC ATA GAC TAG
    assert_equal 'SR', @subseq.aa_seq_in_contig_frame(2) # G T[AC AGT CGG N]NN TCA TAG ACT AG
    assert_equal 'TVX', @subseq.aa_seq_in_contig_frame(3) # GT [ACA GTC GGN] NNT CAT AGA CTA G
    assert_equal 'RL', @subseq.aa_seq_in_contig_frame(4) # CTA GTC TAT GAN N[NC CGA CTG T]AC
    assert_equal 'XDC', @subseq.aa_seq_in_contig_frame(5) # C TAG TCT ATG ANN [NCC GAC TGT] AC
    assert_equal 'PT', @subseq.aa_seq_in_contig_frame(6) # CT AGT CTA TGA NN[N CCG ACT GT]A C
  end

  def test_get_na_coord_for_aa_in_contig_frame
    seq = 'GTACAGTCGGNNNTCATAGACTAGGAGCCTGTACA'
    build_dataset(seq)

    build_subseq(11, 28)
    # (2/5 frame)

    # ----------NNNTCATAGACTAGGAGC-------
    # GTACAGTCGGNNNTCATAGACTAGGAGCCTGTACA

    # -------GCTCCTAGTCTATGANNN----------
    # TGTACAGGCTCCTAGTCTATGANNNCCGACTGTAC

    # 1: X S ^ T R S
    # 2: X H R L G
    # 3: X I D ^ E
    # 4: A P S L ^ X
    # 5: L L V Y X X
    # 6: S ^ S M X X

    assert_nil @subseq.get_na_coord_for_aa_in_contig_frame('A', 2)
    assert_equal 25, @subseq.get_na_coord_for_aa_in_contig_frame('E', 1)
    assert_equal 17, @subseq.get_na_coord_for_aa_in_contig_frame('^', 2)
    assert_equal 15, @subseq.get_na_coord_for_aa_in_contig_frame('H', 3)

    assert_equal 17, @subseq.get_na_coord_for_aa_in_contig_frame('M', 4)
    assert_equal 25, @subseq.get_na_coord_for_aa_in_contig_frame('P', 5)
    assert_equal 15, @subseq.get_na_coord_for_aa_in_contig_frame('X', 6)
  end

  protected

  def build_dataset(seq)
    @dataset = TestDataset.new(seq)
    @contig = @dataset.contig
  end

  def build_subseq(left_id, right_id)
    raise 'No contig' unless contig

    @subseq = ContigSubsequence.new(contig, left_id, right_id)
  end
end
