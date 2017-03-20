class HitsPreparer
  attr_accessor :blast_reader

  def initialize(hits_path:, target:, mode:)
    @blast_reader = BlastReader.new hits_path
  end

  def prepare!
    # back translate and save a translated copy
    # group hits by contigs, save each in separate file
    # for each file: sort, cluster and save to separate gff
  end

  protected

  def back_translate

  end

  def split_by_contigs
  end

  def sort_and_cluster
  end
end
