Dir["#{ROOT_PATH}/lib/contig_elements/zoi/*.rb"].each {|file| require file }

module ContigElements
  class Zoi < ContigElement
    class ZoiHitsException < Exception
    end

    include Polycistronic
    include Annotation
    include Gff

    DIRECTIONS = ['+', '-'].freeze
    FRAMES = (1..6).to_a.freeze
    START_CODON = 'M'
    STOP_CODON = '*'

    attr_reader :contig, :validation_errors, :defection_reasons, :extra_data,
      :gene_begin, :gene_end, :source_frame, :raw_gff
    # start - left border
    # finish - right border
    # gene_begin - real begin of annotated gene (by direction)
    # gene_end - real end of annotated gene (by direction)

    def initialize(contig, start, finish, raw_gff)
      @contig = contig
      @start = start
      @finish = finish
      @extra_data = extra_data
      @source_frame = source_frame
      @raw_gff = raw_gff

      @validation_errors = []
      @defection_reasons = []

      @gene_begin = nil
      @gene_end = nil
    end

    def annotated?
      @gene_begin && @gene_end
    end

    def make_invalid!(reason:)
      unless BadTranscriptsLogger.correct_invalidity_reason?(reason)
        raise "Unknown invalidity reason: #{reason}"
      end

      reason = reason.to_sym
      @validation_errors << reason unless @validation_errors.include?(reason)
    end

    def validate
      make_invalid! reason: :short_transcript if short?
      make_invalid! reason: :no_hits if without_hits?
    end

    def valid?
      @validation_errors.empty?
    end

    def invalid?
      !valid?
    end

    def check_defection
      make_defective! reason: :has_no_sl_mappings if sls.empty?
      make_defective! reason: :fused_genes if blast_hits.count > 1
    end

    def defective?
      return false unless valid?
      @defection_reasons.any?
    end

    def make_defective!(reason:)
      unless BadTranscriptsLogger.correct_defection_reason?(reason)
        raise "Unknown defection reason: #{reason}"
      end

      reason = reason.to_sym

      @defective = true
      @defection_reasons << reason unless @defection_reasons.include?(reason)
    end

    def short?
      (finish - start + 1) < Settings.annotator.transcript_min_size
    end

    def without_hits?
      blast_hits.count == 0
    end

    def blast_hits
      @blast_hits ||= begin
        min_intersection_rate = Settings.annotator.zoi_hit_searching_min_intersection_rate

        contig.blast_hits
              .select_intersected([start, finish])
              .select { |h| h.intersection_rate(self) >= min_intersection_rate }
              .sort_by { |h| h.start }
      end
    end

    def blast_hit_begin
      @blast_hit_begin ||= blast_hit_forward? ?
        blast_hits.first.start :
        blast_hits.last.finish
    end

    def blast_hit_end
      @blast_hit_end ||= blast_hit_forward? ?
        blast_hits.last.finish :
        blast_hits.first.start
    end

    ## hits by coordinates

    def left_blast_hit
      blast_hits.first
    end

    def right_blast_hit
      blast_hits.last
    end

    ## hits by direction

    def begin_blast_hit
      blast_hit_forward? ? left_blast_hit : right_blast_hit
    end

    def end_blast_hit
      blast_hit_forward? ? right_blast_hit : left_blast_hit
    end

    def blast_hit_frame
      left_blast_hit.frame
    end

    def blast_hit_forward?
      left_blast_hit.forward?
    end

    def blast_hit_reverse?
      left_blast_hit.reverse?
    end

    def blast_hit_direction
      bleft_blast_hit.direction
    end

    def merged_gene?
      blast_hits.count > 1
    end

    def left_sls_sorted
      @left_sls_sorted ||= begin
        interval = (start-outer_threshold)..(start+inner_threshold)
        sorted_sls_by_interval interval
      end
    end

    def right_sls_sorted
      @right_sls_sorted ||= begin
        interval = (finish-inner_threshold)..(finish+outer_threshold)
        sorted_sls_by_interval interval
      end
    end

    def sls
      @sls ||= ContigElementCollection.new left_sls_sorted + right_sls_sorted
    end

    def begin_torn?
      torn_by_coord = blast_hit_forward? ? 1 : contig.length
      begin_blast_hit.begin == torn_by_coord
    end

    def end_torn?
      torn_by_coord = forward? ? contig.length : 1
      end_blast_hit.end == torn_by_coord
    end

    def torn_by_coord
      if begin_torn?
        begin_blast_hit.begin
      elsif end_torn?
        end_blast_hit.end
      end
    end

    def normalize!
      @start, @finish = [@start, @finish].sort
    end

    def to_range
      start..finish
    end

    protected

    def inner_threshold
      na_len * Settings.annotator.zoi_sl_searching_inner_multiplier.to_f
    end

    def outer_threshold
      Settings.annotator.zoi_sl_searching_outer_threshold.to_f
    end

    def sorted_sls_by_interval(interval)
      sls = contig.sl_mappings.select_intersected([interval.begin, interval.end])
      ContigElementCollection.new sls.sort_by { |e| [-e.coverage, (e.start-center_coord).abs] }
    end

    def center_coord
      (start+finish)/2
    end
  end
end
