module ContigElements
  class Zoi < ContigElement
    module Annotation
      def annotate
        # Not recommended as it overrides the hard coded table
        table = Bio::CodonTable[1]
        table['tga'] = '^'
        table['tag'] = '^'

        aa_sequence = aa_seq(local_frame)
        idx = aa_sequence.index(START_CODON)

        unless idx
          make_invalid! reason: :cannot_detect_start
          BadTranscriptsLogger.add_to_bin self
          return false
        end

        @gene_start = forward? ? start+idx*3+local_frame-1 : finish-idx*3-local_frame+4

        ######################################
        ######## Stop searching BEGIN ########
        ######################################

        raise "Wrong hit cluster: #{best_blast_hit.data.data}" if best_blast_hit.na_len%3 != 0
        hit_borders = [best_blast_hit.start, best_blast_hit.finish].sort
        hit_end = forward? ? hit_borders[1]+1 : hit_borders[0]-1

        zoi_borders = [start, finish].sort
        zoi_end = forward? ? zoi_borders[1] : zoi_borders[0]

        stop_id = nil

        ## First attempt:
        ## closest STOP to real hit end in the interval from real hit end to the end of zoi

        if (forward? && (hit_end < zoi_end)) || (!forward? && (hit_end > zoi_end))
          first_border, second_border = [hit_end, zoi_end].sort
          stop_id = find_stop(first_border, second_border, hit_end)
          # write_log_to_new_gff('results/annotator/tst.gff', first_border, second_border, direction: '+', contig: nil, extra: {})
        end

        ## Second attempt
        ## closest STOP to real hit end in the interval from the gene start to the end of zoi

        # unless stop_id
        #   first_border = @gene_start
        #   second_border = forward? ? finish : start
        #   first_border, second_border = [first_border, second_border].sort

        #   stop_id = find_stop(first_border, second_border, hit_end)
        # end

        if stop_id.nil?
          make_invalid! reason: :cannot_detect_stop
          BadTranscriptsLogger.add_to_bin self
          return false
        end

        ####################################
        ######## Stop searching END ########
        ####################################

        @gene_finish = stop_id

        if (@gene_finish - @gene_start).abs+1 < Settings.annotator.gene_min_size
          make_invalid! reason: :short_gene
          BadTranscriptsLogger.add_to_bin self
          return false
        end

        puts "Gene #{@raw_gff} (valid: #{@valid}, defective: #{@defection_reason}) has incorrect length [#{@gene_start}, #{@gene_finish}] (#{(@gene_finish-@gene_start).abs+1})" if ((@gene_finish-@gene_start).abs+1)%3 != 0

        true
      end

      def find_stop(left_border, right_border, closest_point)
        stop_aa_region = contig.aa_subseq(left_border, right_border, forward: forward?)

        closest_stop_idx = nil
        closest_dist = nil

        stop_aa_region.to_s.split('').each_with_index do |aa, i|
          if aa == '*'
            current_idx = forward? ? left_border+i*3-1 : right_border-i*3+1
            dist = (current_idx - closest_point).abs

            closest_stop_idx ||= i
            closest_dist ||= dist

            if dist < closest_dist
              closest_dist = dist
              closest_stop_idx = i
            end
          end
        end

        closest_stop_idx = forward? ? left_border+closest_stop_idx*3-1 : right_border-closest_stop_idx*3+1 if closest_stop_idx
        closest_stop_idx
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
        @defective_local_frame ||= begin
          # detection by best blast hit
          if defective_direction == '+'
            (best_blast_hit.start-start)%3+1
          else
            (finish-best_blast_hit.finish)%3+4
          end
        end
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

      def has_no_correct_frames?
        @has_no_correct_frames ||= !can_detect_correct_direction? || correct_local_frame.nil? || !in_bh_cluster_frame?
      end

      def can_detect_correct_direction?
        !sl_mapping.nil?
      end
    end
  end
end
