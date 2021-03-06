require_relative 'helpers/file_splitter.rb'
Dir["#{ROOT_PATH}/lib/preparer/*.rb"].each {|file| require file }

class Preparer
  include Paths

  attr_reader :paths

  def initialize
    @genome_path = SettingsHelper.instance.abs_path_for_setting :genome
    @hits_path = SettingsHelper.instance.abs_path_for_setting :blast_hits
    @transcriptome_path = SettingsHelper.instance.abs_path_for_setting :transcriptome
    @sl_mapping_path = SettingsHelper.instance.abs_path_for_setting :sl_mapping
  end

  def prepare!
    create_tmp_folders

    split_hits_by_contigs
    merge_hits
    back_translate_hits
    clean_hits
    keep_best_nonoverlapped_hits
    save_merged_nonoverlapped_hits

    split_transcripts_by_contigs
    split_sl_mapping_by_contigs

    remove_unnecessary_contigs

    clean_transcriptomes if Settings.annotator.clean_transcriptome

    make_gffs
  end

  def create_tmp_folders
    if Preparer.tmp_abs_pathname.exist?
      Preparer.tmp_abs_pathname.rmtree
    else
      Preparer.tmp_abs_pathname.mkpath
    end

    Preparer.contigs_folder_path.mkpath unless Preparer.contigs_folder_path.exist?

    FastaReader.new(@genome_path).each do |record|
      Preparer.contig_folder_path(record.entry_id).mkpath
    end
  end

  def split_hits_by_contigs
    sort_by_target
    split_csv_by_contigs @hits_path, PATHS[:hits_csv]
  end

  def merge_hits
    each_folder('Merging hits') do |folder|
      hits_path = Preparer.hits_csv_path(folder)
      output_path = Preparer.after_merging_hits_path(folder)
      next unless hits_path.exist?

      reader = BlastReader.new hits_path
      reader.merge_hits! output_path: output_path,
                         target: target,
                         max_distance: 512

      # merged hits showing
      merged_path = Preparer.merged_hits_gff_path(folder)
      File.open(merged_path, 'w') do |f|
        reader.each_hit do |hit|
          next unless hit.merged?
          hit.back_translate_coords! target
          f.puts hit.to_gff(target, show_merged: true)
        end
      end

      reader.close
    end
  end

  def back_translate_hits
    each_folder('Back translating hits') do |folder|
      hits_path = Preparer.after_merging_hits_path(folder)
      output_path = Preparer.back_translated_hits_path(folder)
      next unless hits_path.exist?

      reader = BlastReader.new hits_path
      reader.back_translate! output_path: output_path,
                             target: target
      reader.close
    end
  end

  def clean_hits
    each_folder('Cleaning hits') do |folder|
      hits_path = Preparer.back_translated_hits_path folder
      next unless hits_path.exist?

      reader = BlastReader.new(hits_path)
      reader.cache_hits
      reader.hits.keep_if { |h| h.data[:evalue].to_f < Settings.annotator.max_evalue }
      reader.write_to_file
      reader.close
    end
  end

  def keep_best_nonoverlapped_hits
    each_folder('Selecting nonoverlapped best hits') do |folder|
      input_path = Preparer.back_translated_hits_path folder
      output_path = Preparer.clusters_csv_path folder
      next unless input_path.exist?

      Filterers::BestNonoverlappedBlastHits.new(input_path, output_path: output_path, target: target).perform
      `cp #{Preparer.clusters_csv_path folder} #{Preparer.hits_csv_path folder}`
    end
  end

  def save_merged_nonoverlapped_hits
    each_folder('Saving merged nonoverlapped hits') do |folder|
      input_path = Preparer.clusters_csv_path folder
      output_path = Preparer.merged_nonoverlapped_hits_gff_path folder

      next unless File.file?(input_path)

      merged_hits = {}
      File.open(Preparer.merged_hits_gff_path(folder), 'r').each do |l|
        splitted = l.strip().split("\t")
        notes = splitted[-1].split(';').map{|e| e.split('=')}.to_h
        if notes['Parent']
          merged_hits[notes['Parent']] ||= {}
          merged_hits[notes['Parent']][:elements] ||= []
          merged_hits[notes['Parent']][:elements] << l
        else
          merged_hits[notes['ID']] ||= {}
          merged_hits[notes['ID']][:elements] ||= []
          merged_hits[notes['ID']][:main] = l
          merged_hits[notes['ID']][:start] = splitted[3].to_i
          merged_hits[notes['ID']][:finish] = splitted[4].to_i
        end
      end

      merged_hits = merged_hits.map { |k, v| [k[/.+?(?=_\d{2,}_\d$)/], v] }.to_h

      File.open(output_path, 'w') do |f|
        reader = BlastReader.new input_path
        reader.each_hit do |h|
          id = h.attr_by_target(:id, h.opposite_target(target))
          if merged_hits.include?(id)
            hit_start, hit_finish = [:start, :finish].map { |e| h.attr_by_target(e, target) }.sort
            found_start, found_finish = [:start, :finish].map { |e| merged_hits[id][e] }.sort

            merged_hits[id].each { |l| f.puts l } if hit_start == found_start && hit_finish == found_finish
          end
        end
        reader.close
      end
    end
  end

  def split_transcripts_by_contigs
    split_csv_by_contigs @transcriptome_path, PATHS[:transcripts_csv]
  end

  def split_sl_mapping_by_contigs
    len = `wc -l #{@sl_mapping_path}`.split(/\s+/).first.to_i

    pb = ProgressBar.create title: "Splitting #{@sl_mapping_path}",
                            starting_at: 0,
                            total: len

    FileSplitter.new(@sl_mapping_path, delimiter: "\t", col_index: 0).each_heap(progress_bar: pb) do |seqid, heap|
      next unless Preparer.contig_folder_path(seqid).exist?

      File.open(Preparer.sl_mapping_path(seqid), 'w') do |f|
        heap.each { |line| f.puts line }
      end
    end
  end

  def remove_unnecessary_contigs
    each_folder('Removing unnecessary contigs') do |folder|
      transcripts_path = Preparer.transcripts_csv_path(folder)

      if !transcripts_path.exist? || BlastReader.new(transcripts_path).hits_count.zero?
        Preparer.contig_folder_path(folder).rmtree
      end
    end
  end

  ## By size, by intersection, etc.
  def clean_transcriptomes
    each_folder('Cleaning transcriptomes') do |folder|
      cleaner = TranscriptomeCleaner.new folder, target: target
      result = cleaner.clean

      if result[:good] && result[:good].any?
        writer = BlastWriter.new(result[:good])
        writer.write_hits Preparer.transcripts_csv_path(folder)
      end
    end
  end

  def make_gffs
    each_folder('Making gffs') do |folder|
      csv_to_gff Preparer.hits_csv_path(folder), Preparer.hits_gff_path(folder)
      csv_to_gff Preparer.clusters_csv_path(folder), Preparer.clusters_gff_path(folder)
      csv_to_gff Preparer.transcripts_csv_path(folder), Preparer.transcripts_gff_path(folder)
      csv_to_gff Preparer.clusters_csv_path(folder), Preparer.clusters_extended_gff_path(folder), extend_borders: true
    end
  end

  protected

  def sort_by_target
    blast_reader = BlastReader.new @hits_path
    blast_reader.sort_by! [target_hit_key], ouput_path: Preparer.sorted_by_target_path(:hits)
    @hits_path = Preparer.sorted_by_target_path(:hits)

    blast_reader = BlastReader.new @transcriptome_path
    blast_reader.sort_by! [target_hit_key], ouput_path: Preparer.sorted_by_target_path(:transcripts)
    @transcriptome_path = Preparer.sorted_by_target_path(:transcripts)
  end

  def target
    Settings.annotator.blast_hit_target.to_sym
  end

  def contig_folders
    Dir["#{Preparer.contigs_folder_path.to_s}/*"]
  end

  def each_folder(title)
    folders = contig_folders
    pb = ProgressBar.create(title: title, starting_at: 0, total: folders.count)

    folders.each do |folder|
      yield folder
      pb.increment
    end
  end

  def target_hit_key
    BlastHit::TARGET_KEYS[target][:id]
  end

  def split_csv_by_contigs(csv_path, new_csv_filename)
    blast_reader = BlastReader.new csv_path

    pb = ProgressBar.create title: "Splitting #{Pathname.new(csv_path).basename.to_s}",
                            starting_at: 0,
                            total: blast_reader.hits_count

    splitter = FileSplitters::BlastHits.new(csv_path, delimiter: ',', attr_name: target_hit_key)
    splitter.each_heap(progress_bar: pb) do |seqid, hits|
      next unless Preparer.contig_folder_path(seqid).exist?

      File.open(Preparer.contig_folder_path(seqid, filename: new_csv_filename), 'w') do |f|
        f.puts blast_reader.headers.join(blast_reader.delimiter)
        hits.each { |hit| f.puts hit.to_csv }
      end
    end
  end

  def csv_to_gff(csv_path, gff_path, extend_borders: false)
    return false unless Pathname.new(csv_path).exist?

    File.open(gff_path, 'w') do |f|
      f.puts "##gff-version 3"
      BlastReader.new(csv_path).each_hit do |hit|
        hit.extend_borders!(target) if extend_borders
        f.puts hit.to_gff(target, show_extended: extend_borders)
      end
    end
  end
end
