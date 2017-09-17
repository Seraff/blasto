class BlastHit
  module Merging

    def needs_to_merge_with?(other_hit, target:, max_distance:)
      opposite_target = opposite_target(target)
      opposite_keys = opposite_keys target

      return false if detect_frame(target) != other_hit.detect_frame(target)
      return false if data[:sseqid] != other_hit.data[:sseqid] || data[:qseqid] != other_hit.data[:qseqid]
      return false if distance_to_hit(other_hit, target: target) > max_distance

      first_interval = data[opposite_keys[:start]]..data[opposite_keys[:finish]]
      second_interval = other_hit.data[opposite_keys[:start]]..other_hit.data[opposite_keys[:finish]]
      if IntervalsHelper.intersects?(first_interval, second_interval)
        puts "ACHTUNG! Hit #{data} and #{other_hit.data} shoud be merged. But have intersected query coordinates..."
        return false
      end

      true
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
      opposite_keys = opposite_keys target
      opposite_target = opposite_target target

      borders = (borders(target) + other_hit.borders(target)).sort
      start, finish = borders(target)
      data[keys[:start]], data[keys[:finish]] = start > finish ? [borders.last, borders.first] : [borders.first, borders.last]

      @merging_gaps ||= []
      @merging_gaps << borders[1..2]

      borders = (borders(opposite_target) + other_hit.borders(opposite_target)).sort
      start, finish = borders(opposite_target)
      data[opposite_keys[:start]], data[opposite_keys[:finish]] = start > finish ? [borders.last, borders.first] : [borders.first, borders.last]

      true
    end

    def merged?
      @merging_gaps && @merging_gaps.any?
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
