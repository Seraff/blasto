Dir["#{ROOT_PATH}/lib/contig_elements/zoi/*.rb"].each {|file| require file }

module ContigElements
  class Zoi < ContigElement
    include Polycistronic
    include Annotation
    include Gff

    DIRECTIONS = ['+', '-'].freeze
    FRAMES = (1..6).to_a.freeze
    START_CODON = 'M'
    STOP_CODON = '*'

    attr_reader :contig, :validation_errors, :defection_reasons, :extra_data,
      :gene_start, :gene_finish, :source_frame, :raw_gff
    # start - left border
    # finish - right border
    # gene_start - real start (left or right side, according to frame)
    # gene_finish - real finish

    def initialize(contig, start, finish, raw_gff)
      @contig = contig
      @start = start
      @finish = finish
      @extra_data = extra_data
      @source_frame = source_frame
      @raw_gff = raw_gff

      @validation_errors = []
      @defection_reasons = []

      @gene_start = nil
      @gene_finish = nil
    end

    def annotated?
      @gene_start && @gene_finish
    end

    def validate!

    end

    def make_invalid!(reason: nil)
      unless BadTranscriptsLogger.correct_invalidity_reason?(reason)
        raise "Unknown invalidity reason: #{reason}"
      end

      @valid = false
      @validation_errors << reason if reason
    end

    def make_valid!
      @validation_errors = []
      @valid = true
    end

    def valid?
      if @valid.nil?
        @valid = begin
          @validation_errors = []

          if (finish - start + 1) < Settings.annotator.transcriptome_min_size
            @validation_errors << :short
          elsif blast_hits.count == 0
            @validation_errors << :no_hits
          elsif blast_hits.count > 1
            @validation_errors << :more_than_one_hit
          end

          @validation_errors.empty?
        end
      end

      @valid
    end

    def invalid?
      !valid?
    end

    def defective?
      return false unless valid?

      if @defective.nil?
        @defective = begin
          if sls.empty?
            @defection_reasons << :has_no_sl_mappings
          end

          !@defection_reasons.empty?
        end
      end

      @defective
    end

    def make_defective!(reason:)
      unless BadTranscriptsLogger.correct_defection_reason?(reason)
        raise "Unknown defection reason: #{reason}"
      end

      @defective = true
      @defection_reasons << reason
    end

    def blast_hits
      @blast_hits ||= begin
        contig.blast_hits.select_intersected([start, finish])
      end
    end

    def best_blast_hit
      blast_hits.first
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
      torn_by_coord = forward? ? 1 : contig.length
      best_blast_hit.begin == torn_by_coord
    end

    def end_torn?
      torn_by_coord = forward? ? contig.length : 1
      best_blast_hit.end == torn_by_coord
    end

    def torn_by_coord
      if begin_torn?
        best_blast_hit.begin
      elsif end_torn?
        best_blast_hit.end
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
