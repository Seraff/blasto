class BlastHit
  module Merging
    def needs_to_merge_with?(other_hit, target:, max_distance:)
      opposite_target = opposite_target(target)

      return false if detect_frame(opposite_target) != other_hit.detect_frame(opposite_target)
      return false if data[:sseqid] != other_hit.data[:sseqid] || data[:qseqid] != other_hit.data[:qseqid]

      distance_to_hit(other_hit, target: target) <= max_distance
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

    def to_the_right_of?(other_hit, target:)
      keys = detect_keys target

      left, right = [data[keys[:start]], data[keys[:finish]]].sort
      other_left, other_right = [other_hit.data[keys[:start]], other_hit.data[keys[:finish]]].sort

      if right < other_left
        return false
      elsif left > other_right
        return true
      else
        raise "Incorrect hits location comparison"
      end
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
