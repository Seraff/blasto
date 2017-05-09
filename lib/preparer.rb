class Preparer
  attr_accessor :blast_reader
  CSV_HITS_FILENAME = 'hits.csv'
  GFF_HITS_FILENAME = 'hits.gff'
  CLUSTERS_FILENAME = 'clusters.gff'

  def initialize(hits_path:, target:, mode:)
    @blast_reader = BlastReader.new hits_path
  end

  def prepare!
    create_tmp_folders

    back_translate
    sort_by_target
    split_by_contigs
    sort_and_cluster
  end

  def create_tmp_folders
    if tmp_abs_pathname.exist?
      tmp_abs_pathname.rmtree
    else
      tmp_abs_pathname.mkpath
    end

    hits_by_contigs_folder_pathname.mkpath unless hits_by_contigs_folder_pathname.exist?
  end

  def back_translate
    target = Settings.annotator.blast_hit_target
    pb = ProgressBar.create(title: 'Translating hits', starting_at: 0, total: @blast_reader.hits_count)

    @blast_reader.back_translate! output_path: back_translated_path, target: target, progress_bar: pb
    @blast_reader = BlastReader.new back_translated_path
  end

  def sort_by_target
    @blast_reader.sort_by! target_hit_key, ouput_path: sorted_by_target_path
  end

  def split_by_contigs
    blast_reader.cache_hits

    current_hits_heap = []
    current_seqid = nil

    pb = ProgressBar.create(title: 'Splitting file', starting_at: 0, total: blast_reader.hits_count)

    blast_reader.each_hit do |hit|
      pb.increment

      seqid = hit.data[target_hit_key]

      unless current_seqid
        current_hits_heap << hit
        current_seqid = seqid
        next
      end

      if current_seqid == seqid
        current_hits_heap << hit
      else

        #write heap to files
        contig_folder = hits_by_contigs_folder_pathname + Pathname.new(current_seqid)
        contig_folder.mkpath

        out_csv = File.open(contig_folder.join(CSV_HITS_FILENAME), 'w')
        out_csv.puts blast_reader.headers.join(blast_reader.delimiter)

        out_gff = File.open(contig_folder.join(GFF_HITS_FILENAME), 'w')
        out_gff.puts "##gff-version 3"

        current_hits_heap.each do |bh|
          out_csv.puts bh.to_csv
          out_gff.puts bh.to_gff Settings.annotator.blast_hit_target
        end

        out_csv.close
        out_gff.close

        current_hits_heap = []
        current_hits_heap << hit
        current_seqid = seqid
      end
    end
  end

  def sort_and_cluster
    folders = Dir["#{hits_by_contigs_folder_pathname}/*"]
    pb = ProgressBar.create(title: 'Clustering files', starting_at: 0, total: folders.count)

    folders.each do |folder|
      base_path = Pathname.new(folder)

      hits_path = base_path.join(GFF_HITS_FILENAME)
      clusters_path = base_path.join(CLUSTERS_FILENAME)

      GffClusterizer.new(input: hits_path, output: clusters_path, max_distance: Settings.annotator.clustering_max_distance)
                    .cluster_and_merge

      pb.increment
    end

    true
  end

  protected

  def tmp_abs_pathname
    tmp_dir = Pathname.new(Settings.annotator.tmp_dir)
    tmp_dir.absolute? ? tmp_dir : Pathname.new(ROOT_PATH).join(tmp_dir)
  end

  def sorted_by_target_path
    tmp_abs_pathname + Pathname.new('hits_sorted_by_target.csv')
  end

  def hits_by_contigs_folder_pathname
    tmp_abs_pathname + Pathname.new('hits_by_contigs')
  end

  def back_translated_path
    tmp_abs_pathname + Pathname.new('hits_back_translated.csv')
  end

  def sorted_path
    change_path(back_translated_path, append: :sorted)
  end

  def target_hit_key
    BlastHit::TARGET_KEYS[Settings.annotator.blast_hit_target.to_sym][:id]
  end
end
