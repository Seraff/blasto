require 'bio'
require_relative './blast_hit.rb'

class BlastReader
  attr_accessor :headers, :delimiter, :file

  def initialize(file_name, delimiter=',')
    @hits = []
    file_name = prepare_filename file_name
    @file = File.open(file_name, 'r')
    @delimiter = delimiter

    parse_file
    true
  end

  def hits
    @hits
  end

  def to_gff(filename, for_query = false)
    filename = prepare_filename filename
    f = File.open(filename, 'w')

    f.puts "##gff-version 3"
    hits.each do |h|
      f.puts h.to_gff(for_query)
    end

    f.close
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

  def each_hit
    hits.each { |h| yield h }
  end

  protected

  def parse_file
    @headers = parse_string @file.first

    @file.each do |data|
      data = parse_string data
      raise "Wrong blast file format: node #{data[0]}" if data.count != @headers.count

      data_hash = {}
      @headers.each_with_index do |title, i|
        data_hash[title] = data[i]
      end

      @hits << BlastHit.new(data_hash)
    end

    true
  end

  def parse_string(str)
    str.chomp.split(@delimiter).map(&:strip)
  end

  def prepare_filename(name)
    name = "#{Dir.pwd}/#{name}" if name[0] != '/'
    name
  end
end
