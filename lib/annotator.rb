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

  def each_contig(prefixes: [])
    @genome.each do |fasta_format|
      if prefixes.any?
        next if prefixes.select { |pr| fasta_format.entry_id.start_with?(pr) }.empty?
      end

      yield Contig.new(fasta_format)
    end
  end

  def count
    @count ||= begin
      cnt = @genome.count
      @genome.rewind
      cnt
    end
  end

  def annotate
    raise 'Blast hits must be prepared at first' unless hits_prepared

    pb = ProgressBar.create title: "Annotating",
                            starting_at: 0,
                            total: @genome.count

    @genome.rewind

    BadTranscriptsLogger.remove_old_logs

    contigs_ids = Settings.annotator.contigs_for_annotating || []

    each_contig(prefixes: contigs_ids) do |c|
      c.annotate

      # writing final clusters
      # File.open(Preparer.abs_path_for(ContigElements::BlastHitCluster::FINAL_CLUSTERS_FILENAME), 'w') do |f|
      #   c.blast_hit_clusters.each do |cluster|
      #     f.puts cluster.to_gff
      #   end
      # end

      pb.increment
    end

    BadTranscriptsLogger.gather_full_log
    Contig.gather_full_annotation

    BadTranscriptsLogger.print_reasons_stats
    puts "Annotated: " + `wc -l < #{Preparer.abs_path_for(Contig::ANNOTATION_FILENAME)}`
  end
end
