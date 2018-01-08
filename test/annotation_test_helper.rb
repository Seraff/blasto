module AnnotationTestHelper
  def _detect_borders(string)
    result = []
    %w([ ]).each do |b_type|
      idx = string.index(b_type)
      idx += 1 if idx
      result << idx
    end
    result
  end

  def _parse_entry(line)
    result = {}

    struct, description = line.split('#').map(&:strip)
    data = description ? description.split(',').map(&:strip) : []

    result[:struct] = struct
    return result if data.empty?

    result[:type] = data[0].downcase.to_sym

    left, right = _detect_borders(struct)
    if left && right
      result[:left] = left
      result[:right] = right
    end

    data[1..-1].map(&:downcase).each do |lexem|
      frame_regexp = /frame\s([1-6])/
      coverage_regexp = /coverage\s(\d)/
      if lexem =~ frame_regexp
        result[:frame] = lexem.match(frame_regexp)[1].to_i
      elsif lexem =~ coverage_regexp
        result[:coverage] = lexem.match(coverage_regexp)[1].to_i
      end
    end

    result
  end

  def build_dataset_from_string(str)
    str.split("\n").map(&:strip).each do |line|
        next if line.empty?
        data = _parse_entry(line)

        case data[:type]
        when :genome
          build_dataset data[:struct]
        when :zoi
          dataset.add_zoi [data[:left], data[:right]]
        when :hit
          dataset.add_hit [data[:left], data[:right], data[:frame]]
        when :sl
          dataset.add_sl [data[:left], data[:right], data[:coverage]]
        end
    end
  end

  def assert_annotation(string, direction: '+')
    string = string.split('#')[0].strip
    left, right = _detect_borders(string)
    gene_start = direction == '+' ? left : right
    gene_finish = direction == '+' ? right : left

    assert_equal gene_start, zoi.gene_start
    assert_equal gene_finish, zoi.gene_finish
  end
end
