require 'bio'
require_relative './blast_hit.rb'

class BlastReader
  attr_accessor :headers, :delimiter, :file
  attr_reader :hits

  def initialize(file_name, delimiter=',')
    @hits = []
    file_name = prepare_filename file_name
    @file = File.open(file_name, 'r')
    @delimiter = delimiter

    parse_header
    true
  end

  def each_hit
    if @hits.any?
      @hits.each do |hit|
        yield hit
      end
    else
      @file.each do |data|
        yield parse_hit(data)
      end
    end
  end

  def fetch_hits
    parse_header

    @file.each do |data|
      @hits << parse_hit(data)
    end

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

  def hits_count
    @hits_count ||= `cat #{@file.path} | wc -l`.to_i - 1
  end

  def rewind
    @file.rewind
    parse_header
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

  protected

  def parse_header
    return false if @headers
    @headers = parse_string @file.first
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
