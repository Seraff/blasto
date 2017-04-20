class BlastHit
  module Merging

    def needs_to_merge_with?(other_hit, target:, max_distance:)
      opposite_target = opposite_target(target)

      return false if detect_frame(target) != other_hit.detect_frame(target)
      return false if data[:sseqid] != other_hit.data[:sseqid] || data[:qseqid] != other_hit.data[:qseqid]

      distance_to_hit(other_hit, target: target) <= max_distance
    end

    # we are not using this criterion
    def similar_in_size?(other_hit, target:, threshold:)
      opposite_target = opposite_target(target)

      d1 = distance_to_hit other_hit, target: opposite_target
      d2 = distance_to_hit other_hit, target: target

      [d1, d2].min.to_f/[d1, d2].max.to_f <= threshold
     end

    def merge_with(other_hit, target:)
      keys = detect_keys target

      borders = (borders(target) + other_hit.borders(target)).sort

      start, finish = borders(target)

      if start > finish
        data[keys[:start]], data[keys[:finish]] = borders.last, borders.first
      else
        data[keys[:start]], data[keys[:finish]] = borders.first, borders.last
      end

      true
    end

    def distance_to_hit(other_hit, target:)
      borders = (borders(target) + other_hit.borders(target)).sort
      borders[2] - borders[1]
    end

    def borders(target)
      keys = detect_keys target
      [data[keys[:start]], data[keys[:finish]]]
    end
  end
end
