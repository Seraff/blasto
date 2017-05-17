require_relative './blast_hit.rb'

class BlastWriter
  attr_reader :headers, :delimiter, :hits

  def initialize(hits, delimiter = ',')
    @delimiter = delimiter
    @headers = []
    @hits = hits

    if hits.any?
      @headers = hits.first.headers
    end
  end

  def write_hits(path)
    File.open(path, 'w') do |f|
      f.puts @headers.join(@delimiter)

      hits.each{ |h| f.puts h.to_csv(delimiter: delimiter) }
    end
  end
end
