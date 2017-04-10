require 'bio'
require_relative './blast_hit.rb'

class BlastReader
  attr_accessor :headers, :delimiter, :file
  attr_reader :hits, :hits_cached

  def initialize(file_name, delimiter=',')
    @hits = []
    file_name = prepare_filename file_name
    @file = File.open(file_name, 'r')
    @delimiter = delimiter
    @hits_cached = false

    parse_header
    true
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
      parse_header

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
    @hits_count ||= `cat #{@file.path} | wc -l`.to_i - 1
  end

  def rewind
    @file.rewind
    parse_header || @file.readline
    true
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

  ## back translate everything
  def back_translate(output_file, target:, mode:, progress_bar: nil)
    raise "Invalid target" unless %w(query subject).include?(target.to_s)
    raise "Invalid mode" unless %w(genome transcriptome).include?(mode.to_s)

    with_file_rewinded do
      each_hit do |hit|
        hit.back_translate_coords! target
        gff = hit.to_gff target, extra_data_keys: [:evalue]

        # additional gff processing
        gff_array = gff.split("\t")

        case mode.to_sym
        when :genome
          gff_array[0].gsub!(/_\d+\z/, '')
        when :transcriptome
          gff_array[0].gsub!(/_length_.+/, '')
        end

        output_file.puts gff_array.join("\t")

        progress_bar.increment if progress_bar
      end
    end
  end

  ## merge hits from one reference contig
  def merge_hits(output_file, target:, max_distance: 64)
    # create tmp file with sorted hits by target seqid column

  end

  protected

  def parse_header
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
