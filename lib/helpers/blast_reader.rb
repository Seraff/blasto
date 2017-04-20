require 'bio'
require_relative './blast_hit.rb'

class BlastReader
  attr_accessor :headers, :delimiter, :file
  attr_reader :hits, :hits_cached

  def initialize(file_name, delimiter=',')
    file_name = prepare_filename file_name
    @file = File.open file_name, 'r'
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
    @hits_count ||= `cat #{@file.path} | wc -l`.to_i - 1
  end

  def rewind
    @file.rewind
    parse_headers || @file.readline
    true
  end

  def reopen
    @file.close
    @file = File.open @file.path, 'r'
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

  ## back translate everything
  def back_translate(output_file, target:, mode:, progress_bar: nil, extend_borders: false)
    raise "Invalid target" unless %w(query subject).include?(target.to_s)
    raise "Invalid mode" unless %w(genome transcriptome).include?(mode.to_s)

    with_file_rewinded do
      each_hit do |hit|
        hit.extend_borders! target if extend_borders
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
  def merge_hits(output_file, target:, max_distance: 128, progress_bar: nil)
    prefix = target.to_s[0] # 'q' or 's'

    sort_by! qseqid: :string, sseqid: :string, "#{prefix}frame" => :digit, "#{prefix}start" => :digit

    merged_hits = []
    merging_group = []

    hits.each do |hit|
      if merging_group.empty?
        merging_group << hit

      elsif hit.needs_to_merge_with? merging_group.last, target: target, max_distance: max_distance
        merging_group << hit

      else
        # merge all hits in group
        if merging_group.count == 1
          merged_hits << merging_group.first
        else
          merged_hits << merge_hit_group(merging_group, target: target)
        end
        merging_group = [hit]
      end
      progress_bar.increment if progress_bar
    end

    output_file.puts headers.join(delimiter)

    merged_hits.each do |hit|
      output_file.puts hit.to_blast_row(delimiter: delimiter)
    end

    output_file.close
  end

  def sort_by!(keys)
    indexes = keys.map { |k, type| [@headers.index(k.to_s), type] }.to_h
    raise 'Incorrect keys' if indexes.any?(&:nil?)

    indexes = indexes.map { |i, type| [i+1, type] }.to_h
    index_str = indexes.map { |i, type| "-k #{i},#{i}#{ type == :digit ? 'n' : '' }" }.join ' '
    tmp_path = change_path(file.path, append: 'tmp')

    `(head -n 1 #{file.path} && (tail -n +2 #{file.path} | sort -t#{delimiter} #{index_str})) > #{tmp_path}`
    `mv #{tmp_path} #{file.path}`

    reopen
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

  def merge_hit_group(hits, target:)
    puts "Merging hits group: query #{hits.first.data[:qseqid]}, subject #{hits.first.data[:sseqid]}"
    hit = hits.first
    hits = hits[1..-1]

    hits.each do |h|
      hit.merge_with h, target: target
    end

    hit
  end
end
