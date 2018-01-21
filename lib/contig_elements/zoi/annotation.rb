module ContigElements
  class Zoi < ContigElement
    module Annotation
      HIT_BORDERS_THRESHOLD = 10
      HIT_END_STOP_THRESHOLD = 3

      def annotate
        return false unless valid?

        if begin_torn?
          @gene_begin = torn_by_coord
          make_defective! reason: :torn
        else
          @gene_begin = detect_gene_begin
          return false unless @gene_begin
        end

        if end_torn?
          @gene_end = torn_by_coord
          make_defective! reason: :torn
        else
          @gene_end = detect_gene_end
          return false unless @gene_end
        end

        if (@gene_end - @gene_begin).abs < Settings.annotator.gene_min_size
          make_invalid! reason: :short_gene
        end

        true
      end

      def detect_gene_begin
        if correct_sl
          subs = subsequence_for_sl(correct_sl)
          result = subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)

          raise "Cannot detect start, error in program! #{self.inspect}" unless result
        else
          interval = forward? ?
            [self.begin, blast_hit_begin+HIT_BORDERS_THRESHOLD] :
            [blast_hit_begin-HIT_BORDERS_THRESHOLD, self.begin]

          subs = contig.subsequence(*interval.sort)
          result = subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)

          unless result
            make_invalid! reason: :cannot_detect_start
            return
          end

          make_defective! reason: :hit_in_another_frame
        end

        result
      end

      def detect_gene_end
        _hit_end = forward? ?
          blast_hit_end - HIT_END_STOP_THRESHOLD :
          blast_hit_end + HIT_END_STOP_THRESHOLD
        _end = self.end

        result = nil

        if (forward? && _hit_end < _end) || (reverse? && _hit_end > _end)
          subs = contig.subsequence(*[_hit_end, _end].sort)
          result = subs.get_na_coord_for_aa_in_contig_frame(STOP_CODON, frame)

          if result
            result -= 1 if forward?
            result += 1 if reverse?
          end
        end

        unless result
          make_invalid! reason: :cannot_detect_stop
          return false
        end

        result
      end

      def correct_sl
        @correct_sl ||= begin
          sls = forward? ? left_sls_sorted : right_sls_sorted

          sls.detect do |sl|
            subs = subsequence_for_sl(sl)
            next unless subs
            subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)
          end
        end
      end

      def subsequence_for_sl(sl)
        sl_center = sl.center_coord_for_frame(frame)

        interval = forward? ?
          [sl_center, blast_hit_begin+HIT_BORDERS_THRESHOLD] :
          [blast_hit_begin-HIT_BORDERS_THRESHOLD, sl_center]

        return if interval[0] >= interval[1]

        contig.subsequence(*interval)
      end

      def hit
        left_blast_hit
      end

      def frame
        hit.frame
      end

      def forward?
        blast_hit_forward?
      end

      def reverse?
        blast_hit_reverse?
      end

      def direction
        forward? ? '+' : '-'
      end

      def begin
        forward? ? start : finish
      end

      def end
      forward? ? finish : start
    end
  end
end
end
