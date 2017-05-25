class BadTranscriptsLogger
	LOG_FILENAME = 'bad_transcripts.gff'

	class << self
		def add_to_bin(contig_name, raw_gff: nil, hit: nil, reason: )
			gff = raw_gff || hit.to_gff(Settings.annotator.blast_hit_target.to_sym)
			raise 'Cannot log an empty gff' if gff.nil?

			gff += "_(#{reason})"
			gff += ';color=#CC070E'

			log_file = Preparer.contig_folder_path(contig_name, filename: LOG_FILENAME)
			`touch #{log_file}` unless log_file.exist?

			File.open(log_file, 'a') do |f|
				f.puts gff
			end
		end

		def add_hit_collection_to_bin(contig_name, hits:, reason:)
			hits.each { |h| BadTranscriptsLogger.add_to_bin(contig_name, hit: h, reason: reason) }
		end

		def gather_full_log
			full_log_path = Preparer.abs_path_for(LOG_FILENAME)
      full_log_path.rmtree if full_log_path.exist?

      `touch #{full_log_path}`

      Dir["#{Preparer.contigs_folder_path.to_s}/*"].each do |folder_name|
        path = Preparer.contig_folder_path(folder_name) + Pathname.new(LOG_FILENAME)
        next unless path.exist?

        `cat #{path} >> #{full_log_path}`
      end
		end
	end
end
