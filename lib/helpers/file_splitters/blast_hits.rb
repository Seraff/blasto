require_relative '../blast_hit.rb'

module FileSplitters
  class BlastHits < FileSplitter
    def initialize(file_path, delimiter:, attr_name:)
      File.open(file_path, 'r') { |f| @headers = f.readline.strip.split(delimiter) }
      super file_path, delimiter: delimiter, col_index: @headers.index(attr_name.to_s)
      skip_lines(1)
    end

    def each_heap(progress_bar: nil)
      super do |val, heap|
        hits = heap.map do |raw_hit|
          splitted_hit = raw_hit.split(@delimiter)
          data = @headers.map { |h| [h, splitted_hit[@headers.index(h)]] }.to_h
          ::BlastHit.new(@headers, data)
        end

        yield val, hits
      end
    end
  end
end
