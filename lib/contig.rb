require_relative 'contig_element_collection'

class Contig
  attr_reader :blast_hits, :blast_hit_clusters, :transcripts, :sl_mappings

  ANNOTATION_FILENAME = 'annotation.gff'

  def initialize(fasta_format)
    @fasta = fasta_format

    # caching
    blast_hits
    blast_hit_clusters
    sl_mappings
    zoi
  end

  def seq
    @fasta.seq
  end

  def aa_seq(frame)
    seq.translate frame
  end

  def title
    @fasta.entry_id
  end

  def length
    @fasta.length
  end

  def annotate
    puts "annotating contig #{title}"

    File.open(Preparer.contig_folder_path(title, filename: ANNOTATION_FILENAME), 'w') do |f|
      zoi.each do |z|
        puts "valid? #{z.valid?}"
        f.puts z.to_gff if z.valid?
      end
    end

    true
  end

  def inspect
    denied_vars = [:@blast_hits, :@blast_hit_clusters, :@transcripts, :@sl_mappings, :@fasta]
    vars = instance_variables.select { |v| !denied_vars.include? v }
                             .map { |v| [v, instance_variable_get(v)] }
                             .to_h
                             .map { |k, v| "#{k}='#{v}'"}
                             .join ' '
    "<##{self.class}:#{object_id.to_s(16)} #{vars}>"
  end


  def blast_hits
    @blast_hits ||= begin
      puts 'blast_hits'
      reader = BlastReader.new(Preparer::hits_csv_path(title))
      reader.cache_hits

      elements = reader.hits.map do |h|
        ContigElement.new h.start(target), h.finish(target), h
      end

      ContigElementCollection.new elements
    end
  end

  def blast_hit_clusters
    @blast_hit_clusters ||= begin
      puts 'blast_hit_clusters'
      elements = []
      File.open(Preparer::hit_clusters_path(title), 'r').each do |line|
        next if line.start_with? '#'

        splitted = line.split "\t"
        start = splitted[3].to_i
        finish = splitted[4].to_i
        frame = splitted[7].to_i

        hits = blast_hits.select_intersected([start, finish])
        hits.keep_if { |h| h.data.detect_frame(target) == frame }

        elements << ContigElement.new(start, finish, hits)
      end

      ContigElementCollection.new elements
    end
  end

  def sl_mappings
    @sl_mappings ||= begin
      elements = []

      File.open(Preparer::sl_mapping_clusters_path(title), 'r').each do |line|
        next if line.start_with? '#'

        splitted = line.split "\t"
        start = splitted[1].to_i
        finish = splitted[2].to_i
        elements << ContigElement.new(start, finish, nil)
      end

      ContigElementCollection.new elements
    end
  end

  def zoi
    @zoi ||= begin
      puts 'preparing zoi'
      elements = []

      File.open(Preparer::transcripts_gff_path(title), 'r').each do |line|
        next if line.start_with? '#'

        splitted = line.split "\t"
        start = splitted[3].to_i
        finish = splitted[4].to_i

        elements << Zoi.new(self, start, finish)
      end

      ContigElementCollection.new elements
    end
  end

  def target
    Settings.annotator.blast_hit_target.to_sym
  end

  def detect_orf_frame(na_start, na_stop, direction)
    frames = case direction
    when '+'
      [1, 2, 3]
    when '-'
      [4, 5, 6]
    end

    indexes = frames.map { |f| aa_seq(f).index('M') }
    min_start_id = indexes.index(indexes.compact.min)
    frames[min_start_id]
  end
end
