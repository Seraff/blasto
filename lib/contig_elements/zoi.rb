Dir["#{ROOT_PATH}/lib/contig_elements/zoi/*.rb"].each {|file| require file }

module ContigElements
	class Zoi < ContigElement
		include Polycistronic
		include Annotation

		DIRECTIONS = ['+', '-'].freeze
		FRAMES = (1..6).to_a.freeze
		START_CODON = 'M'

		VALID_COLOR = '#009900'
		DEFECTIVE_COLOR = '#d1b204'
		INVALID_COLOR = '#CC070E'

		attr_reader :contig, :frame,
								:validation_error, :extra_data,
								:gene_start, :gene_finish, :source_frame, :raw_gff

		# start - left border
		# finish - right border
		# gene_start - real start (left or right side, according to frame)
		# gene_finish - real finish

		def initialize(contig, start, finish, raw_gff)
			@contig = contig
			@start = start
			@finish = finish
			@extra_data = extra_data
			@source_frame = source_frame
			@raw_gff = raw_gff

			@validation_error = nil
			@defection_reason = nil

			@gene_start = nil
			@gene_finish = nil
		end

		def annotated?
			!@gene_start.nil? && !@gene_finish.nil?
		end

		def make_invalid!(reason: nil)
			@valid = false
			@validation_error = reason if reason
		end

		def make_valid!
			@valid = true
		end

		def valid?
			if @valid.nil?
				@valid = begin
					@validation_error = nil

					if hit_clusters.count > 1
						@validation_error = :hit_clusters_more_than_one
					elsif sl_mapping.nil?
						@defection_reason ||= :has_no_sl_mappings
					elsif correct_local_frame.nil?
						@defection_reason ||= :cannot_detect_frame
					elsif !in_bh_cluster_frame?
						@defection_reason ||= :hit_cluster_has_another_frame
					end

					@validation_error.nil?
				end
			end

			@valid
		end

		def invalid?
			!valid?
		end

		def defective?
			valid? && !@defection_reason.nil?
		end

		def make_defective!(reason:)
			@defection_reason = reason
		end

		def hit_clusters
			@hit_clusters ||= begin
				contig.blast_hit_clusters.select_intersected([start, finish])
			end
		end

		def hit_cluster
			hit_clusters.first
		end

		def blast_hits
			@blast_hits ||= begin
				contig.blast_hits.select_intersected([start, finish])
			end
		end

		def best_blast_hit
			hit_cluster && hit_cluster.best_blast_hit
		end

		def sl_mappings
			@sl_mappings ||= begin
				zones = [[start-outer_threshold, start+inner_threshold],
								 [finish-inner_threshold, finish+outer_threshold]]
				contig.sl_mappings.select_intersected(*zones)
			end
		end

		def all_sl_mappings
			@inner_sl_mappings ||= begin
				zone = [start-outer_threshold, finish+outer_threshold]
				contig.sl_mappings.select_intersected(zone)
			end
		end

		def sl_mapping
			sl_mappings.first
		end

		def to_gff
			color = if valid?
								if defective?
									DEFECTIVE_COLOR
								else
									VALID_COLOR
								end
							else
							 INVALID_COLOR
							end


			if annotated?
				left, right = [gene_start, gene_finish].sort
				f = global_frame
				d = direction
				notes = "ID=#{SecureRandom.hex}"
				notes += "_#{@defection_reason}" if defective?
			else
				left, right = [start, finish].sort
				f = raw_gff_hash[:frame]
				d = raw_gff_hash[:direction]
				notes = raw_gff_hash[:notes]
			end

			unless valid?
				notes.gsub!(/ID=(.+?)(?=(;|\z))/, "ID=\\1_(#{validation_error});")
			end

			notes += ";color=#{color}"

			if valid?
				if d == '+'
					bbh_finish = best_blast_hit.finish
					extended_bbh_finish = best_blast_hit.extended_finish
				else
					bbh_finish = best_blast_hit.start
					extended_bbh_finish = best_blast_hit.extended_start
				end

				stop_distance = (gene_finish - bbh_finish).abs
				extended_stop_distance = (gene_finish - extended_bbh_finish).abs

				notes += ";distance_to_hit_finish=#{stop_distance}"
				notes += ";distance_to_hit_extended_finish=#{extended_stop_distance}"
				notes += ";bbh_evalue=#{best_blast_hit.data.data[:evalue]}"
				notes += ";bbh_name=#{best_blast_hit.data.data[:qseqid]}"
			end

			f = [4,5,6].include?(f) ? (f - 4) : (f - 1)
			[contig.title, :blast, :gene, left, right, '.', d, f, notes].join("\t")
		end

		def to_gff_as_is
			left, right = [start, finish].sort
			[contig.title, :blast, :gene, left, right, '.', '+', '1', 'yay'].join("\t")
		end

		def normalize!
			@start, @finish = [@start, @finish].sort
		end

		def to_range
			start..finish
		end

		def id
			raw_gff_hash[:notes][/ID=(.+?)(?=(;|\z))/, 1]
		end

		protected

		def raw_gff_hash
			@raw_gff_hash ||= begin
				s = raw_gff.split("\t")
				{ direction: s[6], frame: s[7].to_i, notes: s[8] }
			end
		end

		def in_bh_cluster_frame?
			return false unless hit_cluster

			hit_cluster.extra_data[:frame].to_i == correct_global_frame
		end

		def inner_threshold
			Settings.annotator.zoi_sl_searching_inner_threshold.to_i
		end

		def outer_threshold
			Settings.annotator.zoi_sl_searching_outer_threshold.to_i
		end
	end
end
