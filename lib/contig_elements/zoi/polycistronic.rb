module ContigElements
  class Zoi < ContigElement
    module Polycistronic
      attr_accessor :polycistronic_subzois, :polycistronic_parent

      def polycistronic?
        polycistronic_subzois.any?
      end

      def produced_by_splitting?
        !polycistronic_parent.nil?
      end

      def polycistronic_subzois
        @polycistronic_subzois ||= detect_polycistronic_subzois
      end

      def detect_polycistronic_subzois
        return [] if polycistronic_cutting_places.empty?

        new_zois = []
        last_coord = start

        begin
          polycistronic_cutting_places.each do |coord|
            new_start = last_coord == start ? start : last_coord+1
            new_zois << get_polycistronic_subzoi(new_start, coord)

            last_coord = coord
          end

          new_zois << get_polycistronic_subzoi(last_coord+1, finish) if last_coord < finish
        end

        new_zois
      end

      # cutting places - array of coordinates to cut
      def polycistronic_cutting_places
        @polycistronic_cutting_places ||= begin
          return [] if blast_hits.count <= 1

          groups = blast_hit_groups
          return [] if groups.count <= 1

          # ---[***]------[***]---
          # ----------SL----------

          cutting_places = []

          groups.each_with_index do |group, i|
            next_group = groups[i+1]
            break if next_group.nil?

            left = group.max { |e| e.finish }.finish# - Settings.annotator.polycistronic_sl_threshold
            right = next_group.min { |e| e.start }.start# + Settings.annotator.polycistronic_sl_threshold

            sl = contig.sl_mappings.dup.select_intersected([left, right]).first

            if sl
              cutting_places << sl.start
            else
              cutting_places << (left+right)/2
            end
          end

          cutting_places
        end
      end

      def blast_hit_groups
        groups = []
        sorted = blast_hits.dup.sort_by { |h| h.start }

        current_group = []

        sorted.each do |e|
          found = false

          if current_group.empty?
            current_group << e
            next
          end

          last_in_group = current_group.last
          l = last_in_group.finish
          r = e.start

          sls = contig.sl_mappings.dup.select_intersected([l, r])
          sl_min_coverage = Settings.annotator.sl_polycistrony_min_coverage
          sls = sls.select { |sl| sl.coverage > sl_min_coverage }

          if e.frame == last_in_group.frame && sls.empty? && decent_coverage_in_between?(last_in_group, e)
            # fusion
            current_group << e
            found = true
          end

          unless found
            groups << current_group
            current_group = [e]
          end
        end

        groups << current_group

        groups
      end

      def get_polycistronic_subzoi(start_coord, finish_coord)
        raise "Wrong coords (#{start_coord}, #{finish_coord})" if start_coord > finish_coord

        start_coord = [start, start_coord].sort.last
        finish_coord = [finish, finish_coord].sort.first

        obj = self.class.new(contig, start_coord, finish_coord, raw_gff)
        obj.polycistronic_parent = self
        obj
      end

      def decent_coverage_in_between?(left_el, right_el)
        sorted = [left_el, right_el].sort_by { |h| h.start }

        left_cov = contig.mean_coverage(sorted[0].start, sorted[0].finish)
        right_cov = contig.mean_coverage(sorted[1].start, sorted[1].finish)
        middle_cov = contig.mean_coverage(sorted[0].finish, sorted[1].start)

        min_gene_cov = [left_cov, right_cov].min
        threshold = 0.5 # TODO: to settings?

        middle_cov >= min_gene_cov*threshold
      end
    end
  end
end
