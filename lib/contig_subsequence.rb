class ContigSubsequence
  attr_accessor :left_id, :right_id, :contig

  def initialize(contig, left_id, right_id)
    init_custom_codon_table

    @contig = contig
    @left_id = left_id
    @right_id = right_id
  end

  def seq
    Bio::Sequence::NA.new contig.seq[(left_id-1)..(right_id-1)]
  end

  def get_na_coord_for_aa_in_contig_frame(aa, contig_frame)
    aa_id = aa_seq_in_contig_frame(contig_frame).index(aa)
    return unless aa_id

    subseq_frame = subseq_frame_by_contig_frame(contig_frame)

    if forward_frame?(contig_frame)
      left_id+(subseq_frame-1)+aa_id*3
    else
      right_id-(subseq_frame-4)-aa_id*3
    end
  end

  def aa_seq_in_contig_frame(contig_frame)
    aa_seq_in_frame subseq_frame_by_contig_frame(contig_frame)
  end

  def subseq_frame_by_contig_frame(contig_frame)
    subseq_frame = forward_frame?(contig_frame) ?
      forward_subseq_frame :
      reverse_subseq_frame

    frame_diff = (subseq_frame-contig_frame).abs

    if forward_frame?(contig_frame)
      return 1 if subseq_frame == contig_frame
      contig_frame > subseq_frame ? frame_diff+1 : 3-frame_diff+1
    else
      return 4 if subseq_frame == contig_frame
      contig_frame > subseq_frame ? 4+frame_diff : 4+3-frame_diff
    end
  end

  def aa_seq_in_frame(frame)
    return '' if seq.size < 3
    seq.translate frame
  end

  def forward_subseq_frame
    ((left_id-1)%3)+1
  end

  def reverse_subseq_frame
    4+((contig.length-right_id)%3)
  end

  def forward_frame?(frame)
    if [1, 2, 3].include?(frame)
      true
    elsif [4, 5, 6].include?(frame)
      false
    else
      raise "Wrong frame for ContigSubsequence#forward_frame?(): #{frame}"
    end
  end

  def coords_in_contig_frame(contig_frame)
    frame_in_subs = subseq_frame_by_contig_frame(contig_frame)

    if forward_frame?(contig_frame)
      start = left_id + frame_in_subs - 1
      finish = right_id - ((right_id - start + 1) % 3)
    else
      finish = right_id - (frame_in_subs - 4)
      start = left_id + ((finish - left_id + 1)%3)
    end

    [start, finish]
  end

  def na_len
    seq.length
  end

  protected

  def init_custom_codon_table
    table = Bio::CodonTable[1]
    table['tga'] = '^'
    table['tag'] = '^'
  end
end
