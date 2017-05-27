require_relative 'contig'
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

    pb = ProgressBar.create title: "Annotating",
                            starting_at: 0,
                            total: @genome.count

    @genome.rewind

    BadTranscriptsLogger.remove_old_logs

    each_contig do |c|
      c.annotate
      pb.increment
    end

    BadTranscriptsLogger.gather_full_log
    Contig.gather_full_annotation

    BadTranscriptsLogger.print_reasons_stats
  end
end
