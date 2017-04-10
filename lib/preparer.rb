class Preparer
  attr_accessor :blast_reader
  HITS_FILENAME = 'hits.gff'
  CLUSTERS_FILENAME = 'clusters.gff'

  def initialize(hits_path:, target:, mode:)
    @blast_reader = BlastReader.new hits_path
  end

  def prepare!
    create_tmp_folders
    back_translate
    sort
    split_by_contigs
    sort_and_cluster
  end

  def back_translate
    target = Settings.annotator.blast_hit_target
    mode = Settings.annotator.blast_hit_mode
    outfile = File.open(back_translated_path, 'w')
    pb = ProgressBar.create(title: 'Translating hits', starting_at: 0, total: @blast_reader.hits_count)

    @blast_reader.back_translate outfile, target: target, mode: mode, progress_bar: pb
    outfile.close
  end

  def sort
    `bedtools sort -i #{back_translated_path} > #{sorted_path}`
  end

  def split_by_contigs
    current_data = []
    current_contig = nil

    pb = ProgressBar.create(title: 'Splitting file', starting_at: 0, total: `cat #{sorted_path} | wc -l`.to_i)

    File.open(sorted_path, 'r').each do |entry|
      pb.increment
      next if entry[0] == '#'

      data = entry.split("\t")

      if current_data.empty?
        current_data << data
        current_contig = data[0]
        next
      end

      if current_contig == data[0]
        current_data << data
      else

        contig_folder = hits_by_contigs_folder_pathname + Pathname.new(current_contig)
        contig_folder.mkpath

        f = File.open(contig_folder.join(HITS_FILENAME), 'w')
        f.puts "##gff-version 3"

        current_data.each do |array|
          f.puts array.join("\t")
        end

        current_data = []
        current_data << data
        current_contig = data[0]
      end
    end
  end

  def sort_and_cluster
    folders = Dir["#{hits_by_contigs_folder_pathname}/*"]
    pb = ProgressBar.create(title: 'Clustering files', starting_at: 0, total: folders.count)

    folders.each do |folder|
      base_path = Pathname.new(folder)

      hits_path = base_path.join(HITS_FILENAME)
      clusters_path = base_path.join(CLUSTERS_FILENAME)

      GffClusterizer.new(input: hits_path, output: clusters_path, max_distance: Settings.annotator.clustering_max_distance)
                    .cluster_and_merge

      pb.increment
    end

    true
  end

  def drop_tmp_folders
  end

  protected

  def tmp_abs_pathname
    tmp_dir = Pathname.new(Settings.annotator.tmp_dir)
    tmp_dir.absolute? ? tmp_dir : Pathname.new(ROOT_PATH).join(tmp_dir)
  end

  def hits_by_contigs_folder_pathname
    tmp_abs_pathname + Pathname.new('hits_by_contigs')
  end

  def back_translated_path
    file_name = change_path(@blast_reader.file.path, new_dir: '', new_ext: :gff)
    (tmp_abs_pathname + Pathname.new(file_name).basename).to_s
  end

  def sorted_path
    change_path(back_translated_path, append: :sorted)
  end

  def create_tmp_folders
    tmp_abs_pathname.mkpath unless tmp_abs_pathname.exist?
    hits_by_contigs_folder_pathname.mkpath unless hits_by_contigs_folder_pathname.exist?
  end
end
