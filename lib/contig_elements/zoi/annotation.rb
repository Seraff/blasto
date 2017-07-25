module ContigElements
  class Zoi < ContigElement
    module Annotation
      def annotate
        aa_sequence = aa_seq(local_frame)

        idx = aa_sequence.index(START_CODON)
        unless idx
          make_invalid! reason: :cannot_detect_start
          BadTranscriptsLogger.add_to_bin self
          return false
        end
        @gene_start = forward? ? start+idx*3+local_frame-1 : finish-idx*3-local_frame+4

        raise 'Wrong hit cluster' if best_blast_hit.na_len%3 != 0
        bbh_finish = forward? ? best_blast_hit.finish : best_blast_hit.start
        threshold = (best_blast_hit.na_len/(3*4))*3

        if forward?
          first_stop_border = bbh_finish - threshold + 1
        else
          first_stop_border = bbh_finish + threshold - 1
        end

        second_stop_border = forward? ? finish : start

        first_stop_border, second_stop_border = [first_stop_border, second_stop_border].sort

        stop_na_seq = Bio::Sequence::NA.new contig.seq[first_stop_border-1..second_stop_border-1]
        stop_aa_seq = forward? ? stop_na_seq.translate(1) : stop_na_seq.translate(4)

        hit_border = forward? ? first_stop_border : second_stop_border

        closest_stop_idx = nil
        closest_dist = nil

        bbh_finish = forward? ? best_blast_hit.finish : best_blast_hit.start

        stop_aa_seq.to_s.split('').each_with_index do |aa, i|
          if aa == '*'
            global_idx = forward? ? hit_border+i*3-1 : hit_border-i*3+1
            dist = (global_idx - bbh_finish).abs

            closest_stop_idx ||= i
            closest_dist ||= dist

            if dist < closest_dist
              closest_dist = dist
              closest_stop_idx = i
            end
          end
        end

        if closest_stop_idx.nil?
          make_invalid! reason: :cannot_detect_stop
          BadTranscriptsLogger.add_to_bin self
          return false
        end

        @gene_finish = forward? ?
                       first_stop_border+closest_stop_idx*3-1 :
                       second_stop_border-closest_stop_idx*3+1

        true
      end

      def forward?
        direction == '+'
      end

      def reverse?
        !forward?
      end

      def direction
        if has_no_correct_frames?
          defective_direction
        else
          correct_direction
        end
      end

      def correct_direction
        @correct_direction ||= begin
          middle = (sl_mapping.start+sl_mapping.finish)/2.0
          start_dist = (start - middle).abs
          finish_dist = (finish - middle).abs
          start_dist <= finish_dist ? '+' : '-'
        end
      end

      def defective_direction
        @defective_direction ||= [1,2,3].include?(best_blast_hit.frame.to_i) ? '+' : '-'
      end

      def local_frame
        if has_no_correct_frames?
          defective_local_frame
        else
          correct_local_frame
        end
      end

      def correct_local_frame
        @correct_local_frame ||= begin
          if can_detect_correct_direction?
            frames = correct_direction == '+' ? [1, 2, 3] : [4, 5, 6]

            indexes = frames.map { |f| aa_seq(f).index(START_CODON) }
            min = indexes.compact.min

            frames[indexes.index(min)]
          end
        end
      end

      def defective_local_frame
        @defective_local_frame ||= [1,2,3].include?(defective_global_frame) ? 1 : 4
      end

      def global_frame
        @global_frame ||= begin
          if has_no_correct_frames?
            defective_global_frame
          else
            correct_global_frame
          end
        end
      end

      def correct_global_frame
        @correct_global_frame ||= begin
          local_frame_to_global(correct_local_frame, correct_direction) if correct_local_frame
        end
      end

      def defective_global_frame
        @defective_global_frame ||= best_blast_hit.frame
      end

      def local_frame_to_global(local, dir)
        return nil unless local

        if dir == '+'
          left_idx = start
          frame = (start+(local-1)-1)%3+1
        else
          left_idx = finish
          s = left_idx-(local-4)
          frame = 4+((contig.length-s)%3)
        end

        frame
      end

      def local_frame_to_global_old(local, dir)
        return nil unless local

        left_idx = dir == '+' ? start : finish

        frame_mapping = {
          '+' => {
                    0 => { 1 => 3, 2 => 1, 3 => 2 },
                    1 => { 1 => 1, 2 => 2, 3 => 3 },
                    2 => { 1 => 2, 2 => 3, 3 => 1 }
                 },
          '-' => {
                    0 => { 4 => 6, 5 => 4, 6 => 5 },
                    1 => { 4 => 5, 5 => 6, 6 => 4 },
                    2 => { 4 => 4, 5 => 5, 6 => 6 }
                 }
        }
        frame_mapping[dir][left_idx%3][local]
      end

      def has_no_correct_frames?
        @has_no_correct_frames ||= !can_detect_correct_direction? || correct_local_frame.nil? || !in_bh_cluster_frame?
      end

      def can_detect_correct_direction?
        !sl_mapping.nil?
      end
    end
  end
end
