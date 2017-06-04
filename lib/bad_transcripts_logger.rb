class BadTranscriptsLogger
	LOG_FILENAME = 'bad_transcripts.gff'
	TOTAL_STAT_FILENAME = 'stats.txt'

	REASONS = [
							 :short,
							 :intersected,
							 :merged,
							 :totally_covered,
							 :hit_clusters_more_than_one,
							 :has_no_sl_mappings,
							 :cannot_detect_frame,
							 :hit_cluster_has_another_frame,
							 :no_hit_clusters,
							 :cannot_detect_stop
						]

	class << self
		def remove_old_logs
			path = Preparer.abs_path_for(LOG_FILENAME)
			path.rmtree if path.exist?

			Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
				path = Preparer.contig_folder_path(folder_name, filename: LOG_FILENAME)
				path.rmtree if path.exist?
			end
		end

		def add_to_bin(zoi)
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
				f.puts "total: #{total_count}"
			end

			puts `cat #{path}`
		end
	end
end
