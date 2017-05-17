require_relative 'contig'
require_relative 'frame'
require_relative 'frame_collection'
require_relative 'hit'
require_relative 'hit_cluster'
require_relative 'hit_collection'
require_relative 'preparer'

class Annotator
  attr_accessor :genome, :hits_prepared

  def initialize
    @genome_path = SettingsHelper.instance.abs_path_for_setting :genome

    @reads_path = SettingsHelper.instance.abs_path_for_setting :genome_reads
    @hits_prepared = Settings.annotator.skip_preparation

    @genome = Bio::FlatFile.open(Bio::FastaFormat, @genome_path)
  end

  def prepare
    return true if hits_prepared
    preparer = Preparer.new

    return false unless preparer.prepare!
    @hits_prepared = true
  end

  def each_contig
    @genome.each do |fasta_format|
      yield Contig.new(fasta_format)
    end
  end

  def annotate
    raise 'Blast hits must be prepared at first' unless hits_prepared

    each_contig do |c|
      c.annotate
      break
    end
  end
end
