module ContigElements
  class Zoi < ContigElement
    module Annotation
      HIT_BORDERS_THRESHOLD = 10

      def annotate
        if begin_torn?
          @gene_start = torn_by_coord
          make_defective! reason: :torn
        else
          @gene_start = detect_gene_start
          return false unless @gene_start
        end

        if end_torn?
          @gene_finish = torn_by_coord
          make_defective! reason: :torn
        else
          @gene_finish = detect_gene_finish
          return false unless @gene_finish
        end

        true
      end

      def detect_gene_start
        if correct_sl
          subs = subsequence_for_sl(correct_sl)
          result = subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)

          raise "Cannot detect start, error in program! #{self.inspect}" unless result
        else
          interval = forward? ?
            [self.begin, hit.begin+HIT_BORDERS_THRESHOLD] :
            [hit.begin-HIT_BORDERS_THRESHOLD, self.begin]

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

      def detect_gene_finish
        subs = contig.subsequence(*[hit.end, self.end].sort) #TODO: without sorting, but hit.end - threshold, self.end + threshold
        result = subs.get_na_coord_for_aa_in_contig_frame(STOP_CODON, frame)

        if result
          result -= 1 if forward?
          result += 1 if reverse?
        end

        unless result
          make_invalid! reason: :cannot_detect_stop
          return false
        end

        result
      end

      def correct_sl
        @correct_sl ||= begin
          sls = hit.forward? ? left_sls_sorted : right_sls_sorted

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
          [sl_center, hit.begin+HIT_BORDERS_THRESHOLD] :
          [hit.begin-HIT_BORDERS_THRESHOLD, sl_center]

        return if interval[0] >= interval[1]

        contig.subsequence(*interval)
      end

      def hit
        best_blast_hit
      end

      def frame
        hit.frame
      end

      def forward?
        hit.forward?
      end

      def reverse?
        hit.reverse?
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
