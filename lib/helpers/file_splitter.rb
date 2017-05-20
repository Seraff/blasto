class FileSplitter
  attr_reader :delimiter

  def initialize(file_path, delimiter:, col_index:)
    @file = File.open(file_path, 'r')
    @delimiter = delimiter
    @col_index = col_index
  end

  def skip_lines(n)
    n.times { @file.readline }
  end

  def each_heap(progress_bar: nil)
    current_heap = []
    current_value = nil

    @file.each_line do |line|
      progress_bar.increment if progress_bar

      val = line.split(@delimiter)[@col_index]

      unless current_value
        current_heap << line
        current_value = val
        next
      end

      if current_value == val
        current_heap << line
      else

        yield current_value, current_heap.map(&:strip)

        current_heap = []
        current_heap << line
        current_value = val
      end
    end

    yield current_value, current_heap.map(&:strip)
  end
end

require_relative './file_splitters/blast_hits.rb'
