require_relative 'contig_element_collection'

class Contig
  attr_reader :blast_hits, :blast_hit_clusters, :transcripts, :sl_mappings

  ANNOTATION_FILENAME = 'annotation.gff'

  class << self
    def gather_full_annotation
      full_ann_path = Preparer.abs_path_for(ANNOTATION_FILENAME)
      full_ann_path.rmtree if full_ann_path.exist?
      `touch #{full_ann_path}`

      Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
        path = Preparer.contig_folder_path(folder_name) + Pathname.new(ANNOTATION_FILENAME)
        next unless path.exist?

        `cat #{path} >> #{full_ann_path}`
      end
    end
  end

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
    return false unless Preparer.contig_folder_path(title).exist?

    File.open(Preparer.contig_folder_path(title, filename: ANNOTATION_FILENAME), 'w') do |f|
      zoi.each do |z|
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
      elements = []
      path = Preparer::hits_csv_path(title)

      if path.exist?
        reader = BlastReader.new path
        reader.cache_hits

        elements = reader.hits.map do |h|
          ContigElement.new h.start(target), h.finish(target), h
        end
      end

      ContigElementCollection.new elements
    end
  end

  def blast_hit_clusters
    @blast_hit_clusters ||= begin
      elements = []

      path = Preparer::hit_clusters_path(title)

      if path.exist?
        File.open(path, 'r').each do |line|
          next if line.start_with? '#'

          splitted = line.split "\t"
          start = splitted[3].to_i
          finish = splitted[4].to_i
          frame = splitted[7].to_i

          hits = blast_hits.select_intersected([start, finish])
          hits.keep_if { |h| h.data.detect_frame(target) == frame }

          elements << ContigElement.new(start, finish, hits, extra_data: { frame: frame })
        end
      end

      ContigElementCollection.new elements
    end
  end

  def sl_mappings
    @sl_mappings ||= begin
      elements = []
      path = Preparer::sl_mapping_clusters_path(title)

      if path.exist?
        File.open(path, 'r').each do |line|
          next if line.start_with? '#'

          splitted = line.split "\t"
          start = splitted[1].to_i
          finish = splitted[2].to_i
          elements << ContigElement.new(start, finish, nil)
        end
      end

      ContigElementCollection.new elements
    end
  end

  def zoi
    @zoi ||= begin
      elements = []
      path = Preparer::transcripts_gff_path(title)

      if path.exist?
        File.open(path, 'r').each do |line|
          line = line.strip
          next if line.start_with? '#'

          splitted = line.split "\t"
          start = splitted[3].to_i
          finish = splitted[4].to_i

          elements << Zoi.new(self, start, finish, extra_data: { raw_gff: line })
        end
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
