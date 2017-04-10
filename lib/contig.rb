class Contig
  def initialize(fasta_format)
    @fasta = fasta_format
  end

  def seq
    @fasta.seq
  end

  def title
    @fasta.entry_id
  end

  def annotate
    # TODO: main code here
    # prepare hit clusters
    # choose best hits
    # save a result to gff
  end

  protected

  def build_hit_clusters
    ##
  end
end
