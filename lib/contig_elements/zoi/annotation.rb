\
module ContigElements
  class Zoi < ContigElement
    module Annotation
      HIT_BORDERS_THRESHOLD = 10

      def annotate
        if correct_sl
          subs = subsequence_for_sl(correct_sl)
          @gene_start = subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)

          raise "Cannot detect start, error in program! #{self.inspect}" unless @gene_start
        else
          interval = forward? ?
            [self.begin, hit.begin+HIT_BORDERS_THRESHOLD] :
            [hit.end-HIT_BORDERS_THRESHOLD, self.end]

          subs = contig.subsequence *interval.sort
          @gene_start = subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)

          unless @gene_start
            make_invalid! reason: :cannot_detect_start
            return false
          end

          make_defective! reason: :blast_hit_has_another_frame
        end

        subs = contig.subsequence *[hit.end, self.end].sort
        @gene_finish = subs.get_na_coord_for_aa_in_contig_frame(STOP_CODON, frame)

        if @gene_finish
          @gene_finish -= 1 if forward?
          @gene_finish += 1 if reverse?
        end

        unless @gene_finish
          make_invalid! reason: :cannot_detect_stop
          return false
        end

        true
      end

      def correct_sl
        @correct_sl ||= begin
          sls = hit.forward? ? left_sls_sorted : right_sls_sorted

          sls.detect do |sl|
            subs = subsequence_for_sl(sl)
            subs.get_na_coord_for_aa_in_contig_frame(START_CODON, frame)
          end
        end
      end

      def subsequence_for_sl(sl)
        sl_center = sl.center_coord_for_frame(frame)

        interval = forward? ?
          [sl_center, hit.begin+HIT_BORDERS_THRESHOLD] :
          [hit.end-HIT_BORDERS_THRESHOLD, sl_center]

        return if interval[0] >= interval[1]

        contig.subsequence *interval
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
