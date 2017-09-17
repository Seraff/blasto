require 'bio'
require_relative './blast_hit.rb'
require_relative 'blast_reader/utils.rb'

class BlastReader
  include Utils

  attr_accessor :headers, :delimiter, :file
  attr_reader :hits, :hits_cached

  def initialize(file_name, delimiter=',')
    file_name = prepare_filename file_name.to_s
    @file = File.open file_name.to_s, 'r'
    @delimiter = delimiter
    @hits_cached = false

    parse_headers
    true
  end

  def hits
    return @hits if hits_cached?

    cache_hits
    @hits
  end

  def hits_cached?
    @hits_cached
  end

  def each_hit
    if hits_cached?
      @hits.each do |hit|
        yield hit
      end
    else
      @file.each do |data|
        yield parse_hit(data)
      end
    end
  end

  def cache_hits
    @hits = []

    with_file_rewinded do
      parse_headers

      @file.each do |data|
        @hits << parse_hit(data)
      end
    end

    @hits_cached = true

    true
  end

  def modify_hits(output_path)
    outfile = File.open output_path, 'w'
    outfile.puts @headers.join(@delimiter)

    each_hit do |h|
      result = yield h

      outfile.puts stringify_hit(h) if result
    end

    outfile.close
  end

  def with_file_rewinded
    rewind
    yield
    rewind
  end

  def hits_count
    @hits_count ||= begin
      if hits_cached?
        @hits.count
      else
        `cat #{@file.path} | wc -l`.to_i - 1
      end
    end
  end

  def rewind
    return true if hits_cached?
    @file.rewind
    parse_headers || @file.readline
    true
  end

  def reopen(new_file_path=nil)
    @file.close
    @file = File.open (new_file_path || @file.path), 'r'
    @headers = nil
    parse_headers

    cache_hits if hits_cached?

    true
  end

  def close
    @file.close
    @headers = nil
    @hits_cached = false
  end

  ## prevents printing all the hits in console
  def inspect
    denied_vars = [:@hits]
    vars = instance_variables.select { |v| !denied_vars.include? v }
                             .map { |v| [v, instance_variable_get(v)] }
                             .to_h
                             .map { |k, v| "#{k}='#{v}'"}
                             .join ' '
    "<##{self.class}:#{object_id.to_s(16)} #{vars}>"
  end

  def write_to_file(path=nil)
    path = path || @file.path

    File.open(path, 'w') do |f|
      f.puts @headers.join(@delimiter)
      hits.each{ |h| f.puts h.to_csv(delimiter: delimiter) }
    end
  end

  def raw_headers_string
    headers.join(delimiter)
  end

  def hits=(hits)
    @hits_count = nil
    @hits = hits
  end

  protected

  def parse_headers
    return false if @headers
    @headers = parse_string @file.first
    true
  end

  def parse_hit(data)
    raise 'Header is not parsed' unless @headers

    data = parse_string data
    raise "Wrong blast file format: node #{data[0]}" if data.count != @headers.count

    data_hash = {}
    @headers.each_with_index do |title, i|
      data_hash[title] = data[i]
    end

    BlastHit.new(@headers, data_hash)
  end

  def parse_string(str)
    str.chomp.split(@delimiter).map(&:strip)
  end

  def stringify_hit(hit)
    @headers.map { |h| hit.data[h] }.join @delimiter
  end

  def prepare_filename(name)
    name = "#{Dir.pwd}/#{name}" if name[0] != '/'
    name
  end
end
