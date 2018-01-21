class BadTranscriptsLogger
  LOG_FILENAME = 'bad_transcripts.gff'
  TOTAL_STAT_FILENAME = 'stats.txt'

  INVALIDITY_REASONS = [
    :intersected,
    :merged,
    :intersected_by_siblings,
    :covered_by_siblings,
    :cannot_detect_frame,
    :no_hits,
    :cannot_detect_start,
    :cannot_detect_stop,
    :short_transcript,
    :short_gene,
    :splitted
  ]

  DEFECTION_REASONS = [
    :hit_in_another_frame,
    :has_no_sl_mappings,
    :torn,
    :produced_by_splitting,
    :fused_genes
  ]

  class << self
    def remove_old_logs
      return if BadTranscriptsLogger.test_environment?

      path = Preparer.abs_path_for(LOG_FILENAME)
      path.rmtree if path.exist?

      Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
        path = Preparer.contig_folder_path(folder_name, filename: LOG_FILENAME)
        path.rmtree if path.exist?
      end
    end

    def add_to_bin(zoi)
      return if BadTranscriptsLogger.test_environment?

      log_file = Preparer.contig_folder_path(zoi.contig.title, filename: LOG_FILENAME)
      `touch #{log_file}` unless log_file.exist?

      File.open(log_file, 'a') do |f|
        f.puts zoi.to_gff
      end
    end

    def add_hit_collection_to_bin(zois)
      zois.each { |z| BadTranscriptsLogger.add_to_bin(z) }
    end

    def gather_full_log
      return if BadTranscriptsLogger.test_environment?

      full_log_path = Preparer.abs_path_for(LOG_FILENAME)
      full_log_path.rmtree if full_log_path.exist?

      `touch #{full_log_path}`

      Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
        path = Preparer.contig_folder_path(folder_name, filename: LOG_FILENAME)
        next unless path.exist?

        `cat #{path} >> #{full_log_path}`
      end
    end

    def print_reasons_stats
      return if BadTranscriptsLogger.test_environment?

      puts
      path = Preparer.abs_path_for(TOTAL_STAT_FILENAME)
      annotation_path = Preparer.abs_path_for(Contig::ANNOTATION_FILENAME)
      path.rmtree if path.exist?

      File.open(path, 'w') do |f|
        total_count = 0

        stats = { reasons: {}, total: {} }

        { invalid: INVALIDITY_REASONS, defections: DEFECTION_REASONS }.each do |key, vals|
          stats[:reasons][key] = []
          stats[:total][key] = 0

          vals.each do |reason|
            count = `cat #{annotation_path} | grep #{reason} | wc -l`.to_i
            stats[:reasons][key] << "\t#{reason}: #{count}" if count > 0
            stats[:total][key] += count
          end
        end

        stats[:reasons].each do |group, vals|
          next if vals.empty?
          f.puts "#{group}:"
          vals.each { |v| f.puts v }
          f.puts
        end

        valid_count = `grep -c 'color=#{ContigElements::Zoi::Gff::VALID_COLOR}' #{annotation_path}`.to_i
        defective_count = `grep -c 'color=#{ContigElements::Zoi::Gff::DEFECTIVE_COLOR}' #{annotation_path}`.to_i
        invalid_count = `grep -c 'color=#{ContigElements::Zoi::Gff::INVALID_COLOR}' #{annotation_path}`.to_i

        f.puts "invalid: #{invalid_count}"
        f.puts "defective: #{defective_count}"
        f.puts "valid: #{valid_count}"
        f.puts
        f.puts "total good: #{valid_count+defective_count}"
        f.puts "total: #{`cat #{annotation_path} | wc -l`}"
        f.puts
      end

      puts `cat #{path}`
    end

    def test_environment?
      ENV['BLASTO_ENV'] == 'TEST'
    end

    def correct_invalidity_reason?(reason)
      INVALIDITY_REASONS.include? reason
    end

    def correct_defection_reason?(reason)
      DEFECTION_REASONS.include? reason
    end
  end
end
