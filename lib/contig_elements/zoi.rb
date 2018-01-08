Dir["#{ROOT_PATH}/lib/contig_elements/zoi/*.rb"].each {|file| require file }

module ContigElements
  class Zoi < ContigElement
    include Polycistronic
    include Annotation

    DIRECTIONS = ['+', '-'].freeze
    FRAMES = (1..6).to_a.freeze
    START_CODON = 'M'
    STOP_CODON = '*'

    VALID_COLOR = '#009900'
    DEFECTIVE_COLOR = '#d1b204'
    INVALID_COLOR = '#CC070E'

    attr_reader :contig, :validation_error, :defection_reason, :extra_data,
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

      @validation_error = nil
      @defection_reason = nil

      @gene_start = nil
      @gene_finish = nil
    end

    def annotated?
      !@gene_start.nil? && !@gene_finish.nil?
    end

    def validate!

    end

    def make_invalid!(reason: nil)
      @valid = false
      @validation_error = reason if reason
    end

    def make_valid!
      @valid = true
    end

    def valid?
      if @valid.nil?
        @valid = begin
          @validation_error = nil

          if (finish - start + 1) < Settings.annotator.transcriptome_min_size
            @validation_error = :short
          elsif blast_hits.count == 0
            @validation_error = :no_hits
          elsif blast_hits.count > 1
            @validation_error = :more_than_one_hit
          end

          @validation_error.nil?
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
            @defection_reason = :has_no_sl_mappings
          end
        end

        !@defection_reason.nil?
      end
    end

    def make_defective!(reason:)
      @defection_reason = reason
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

    def to_gff
      color = if valid?
        if defective?
          DEFECTIVE_COLOR
        else
          VALID_COLOR
        end
      else
        INVALID_COLOR
      end


      if annotated?
        left, right = [gene_start, gene_finish].sort
        f = frame
        d = direction
        notes = "ID=#{SecureRandom.hex}"
        notes += "_#{@defection_reason}" if defective?
      else
        left, right = [start, finish].sort
        f = raw_gff_hash[:frame]
        d = raw_gff_hash[:direction]
        notes = raw_gff_hash[:notes]
      end

      if invalid?
        invalid_hash = SecureRandom.hex
        notes.gsub!(/ID=(.+?)(?=(;|\z))/, "ID=#{invalid_hash}_(#{validation_error});")
      end

      notes += ";color=#{color}"

      if valid?
        if d == '+'
          bbh_finish = best_blast_hit.finish
          extended_bbh_finish = best_blast_hit.extended_finish
        else
          bbh_finish = best_blast_hit.start
          extended_bbh_finish = best_blast_hit.extended_start
        end

        stop_distance = (gene_finish - bbh_finish).abs
        extended_stop_distance = (gene_finish - extended_bbh_finish).abs

        notes += ";distance_to_hit_finish=#{stop_distance}"
        notes += ";distance_to_hit_extended_finish=#{extended_stop_distance}"
        notes += ";bbh_evalue=#{best_blast_hit.data.data[:evalue]}"
        notes += ";bbh_name=#{best_blast_hit.data.data[:qseqid]}"
      end

      f = [4,5,6].include?(f) ? (f - 4) : (f - 1)
      [contig.title, :blast, :gene, left, right, '.', d, f, notes].join("\t")
    end

    def to_gff_as_is
      left, right = [start, finish].sort
      [contig.title, :blast, :gene, left, right, '.', '+', '1', 'yay'].join("\t")
    end

    def normalize!
      @start, @finish = [@start, @finish].sort
    end

    def to_range
      start..finish
    end

    def id
      raw_gff_hash[:notes][/ID=(.+?)(?=(;|\z))/, 1]
    end

    protected

    def raw_gff_hash
      @raw_gff_hash ||= begin
        s = raw_gff.split("\t")
        { direction: s[6], frame: s[7].to_i, notes: s[8] }
      end
    end

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
