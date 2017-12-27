class BadTranscriptsLogger
	LOG_FILENAME = 'bad_transcripts.gff'
	TOTAL_STAT_FILENAME = 'stats.txt'

	REASONS = [
							 :short_transcript,
							 :intersected,
							 :merged,
							 :totally_covered,
							 :hit_clusters_more_than_one,
							 :has_no_sl_mappings,
							 :cannot_detect_frame,
							 :hit_cluster_has_another_frame,
							 :no_hit_clusters,
							 :cannot_detect_start,
							 :cannot_detect_stop,
							 :short_gene
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

		def gather_all_clusters
			return if BadTranscriptsLogger.test_environment?

			clusters_filename = 'hit_clusters.gff'
			full_clusters_path = Preparer.abs_path_for(clusters_filename)
      full_clusters_path.rmtree if full_clusters_path.exist?

      `touch #{full_clusters_path}`

      Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
        path = Preparer.contig_folder_path(folder_name, filename: clusters_filename)
        next unless path.exist?

        `cat #{path} | grep "^[^#;]" >> #{full_clusters_path}`
      end
		end

		def print_reasons_stats
			return if BadTranscriptsLogger.test_environment?

			puts
			path = Preparer.abs_path_for(TOTAL_STAT_FILENAME)
			path.rmtree if path.exist?

			File.open(path, 'w') do |f|
				total_count = 0

				REASONS.each do |reason|
					count = `cat #{Preparer.abs_path_for(LOG_FILENAME)} | grep #{reason} | wc -l`.to_i
					total_count += count
					f.puts "#{reason}: #{count}"
				end

				f.puts
				f.puts "invalid: #{total_count}"
			end

			puts `cat #{path}`
		end

		def test_environment?
			ENV['BLASTO_ENV'] == 'TEST'
		end
	end
end
