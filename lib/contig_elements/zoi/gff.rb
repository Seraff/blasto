module ContigElements
  class Zoi < ContigElement
    module Gff
      VALID_COLOR = '#009900'
      DEFECTIVE_COLOR = '#d1b204'
      INVALID_COLOR = '#CC070E'

      def to_gff
        if annotated?
          left, right = [gene_begin, gene_end].sort
          f = frame
          d = direction
        else
          left, right = [start, finish].sort
          f = raw_gff_hash[:frame]
          d = raw_gff_hash[:direction]
        end

        notes = "ID=#{SecureRandom.hex}"

        if defective?
          notes += ";defections=#{defection_reasons_as_string}"
        end

        if invalid?
          notes += ";invalidity=#{validation_errors_as_string}"
        end

        notes += ";color=#{color}"

        # hits info
        if valid?
          blast_hits.each_with_index do |hit, i|
            notes += ";hit_#{i}_evalue=#{hit.evalue}"
            notes += ";hit_#{i}_seqid=#{hit.seqid}"
          end

          extended_bbh_finish = forward? ?
            end_blast_hit.extended_finish :
            end_blast_hit.extended_start
          begin
            stop_distance = (gene_end - end_blast_hit.end).abs
          rescue NoMethodError
            binding.pry
          end
          extended_stop_distance = (gene_end - extended_bbh_finish).abs

          notes += ";distance_to_hit_finish=#{stop_distance}"
          notes += ";distance_to_hit_extended_finish=#{extended_stop_distance}"
        end

        status = if valid? || defective? || canot_detect_all_borders?
          :good
        else
          :bad
        end

        notes += ";status=#{status}"

        f = [4,5,6].include?(f) ? (f - 4) : (f - 1)
        [contig.title, :blast, :gene, left, right, '.', d, f, notes].join("\t")
      end

      def to_gff_as_is
        left, right = [start, finish].sort
        [contig.title, :blast, :gene, left, right, '.', '+', '1', 'yay'].join("\t")
      end

      def validation_errors_as_string
        validation_errors.map(&:to_s).join(',')
      end

      def color
        if valid?
          if defective?
            DEFECTIVE_COLOR
          else
            VALID_COLOR
          end
        else
          INVALID_COLOR
        end
      end

      def defection_reasons_as_string
        defection_reasons.map(&:to_s).join(',')
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
    end
  end
end
