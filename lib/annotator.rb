require_relative 'contig'
require_relative 'frame'
require_relative 'frame_collection'
require_relative 'hit'
require_relative 'hit_cluster'
require_relative 'hit_collection'
require_relative 'preparer'

class Annotator
  attr_accessor :genome, :hits_prepared

  def initialize(genome_path:, hits_path:, reads_path:, skip_preparation: false)
    @genome_path = genome_path
    @hits_path = hits_path
    @reads_path = reads_path
    @hits_prepared = skip_preparation

    @genome = Bio::FlatFile.open(Bio::FastaFormat, @genome_path)
  end

  def prepare(target:, mode:)
    return true if hits_prepared
    preparer = Preparer.new hits_path: @hits_path, target: target, mode: mode

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
