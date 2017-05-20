class Contig
  attr_accessor :blast_hits
  # has many Blast Hits
  # has many ZOI
  # has many

  def initialize(fasta_format)
    @fasta = fasta_format

    hits_path = SettingsHelper.instance.tmp_abs_pathname
    hits_path += Pathname.new(Preparer::SPLITTED_DATA_FOLDER)
    hits_path += Pathname.new("#{title}/hits.csv")

    @hits_reader = BlastReader.new hits_path

    init_hits
  end

  def seq
    @fasta.seq
  end

  def title
    @fasta.entry_id
  end

  def annotate
    puts "annotating contig #{title}"
    # make ZOI collection
  end

  protected

  def init_hits

  end
end
