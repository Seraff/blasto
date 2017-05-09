class Contig
  attr_accessor :blast_hits
  # has many Blast Hits
  # has many ZOI
  # has many

  def initialize(fasta_format)
    @fasta = fasta_format
    init_hits
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

  def init_hits

  end
end
