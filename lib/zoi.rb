class Zoi
	DIRECTIONS = ['+', '-'].freeze
	FRAMES = (1..6).to_a.freeze
	START_CODON = 'M'

	VALID_COLOR = '#009900'
	INVALID_COLOR = '#CC070E'

	attr_reader :contig, :start, :finish, :direction, :frame,
							:validation_error, :extra_data,
							:gene_start, :gene_finish, :source_frame, :raw_gff

	def initialize(contig, start, finish, raw_gff)
		@contig = contig
		@start = start
		@finish = finish
		@extra_data = extra_data
		@source_frame = source_frame
		@raw_gff = raw_gff

		@validation_error = nil
		@gene_start = nil
		@gene_finish = nil
	end

	def annotate
		idx = aa_seq(local_frame).index(START_CODON)

		@gene_start = forward? ? start+idx*3+local_frame-1 : finish-idx*3-local_frame+4
		@gene_finish = forward? ? finish : start #TODO

		self
	end

	def annotated?
		!@gene_start.nil? && !@gene_finish.nil?
	end

	def make_invalid!(reason: nil)
		@valid = false
		@validation_error = reason if reason
	end

	def valid?
		if @valid.nil?
			@valid = begin
				@validation_error = nil

				if hit_clusters.count > 1
					@validation_error = :hit_clusters_more_than_one
				elsif sl_mapping.nil?
					@validation_error = :has_no_sl_mappings
				elsif local_frame.nil?
					@validation_error = :cannot_detect_frame
				elsif !in_bh_cluster_frame?
					@validation_error = :hit_cluster_has_another_frame
				end

				@validation_error.nil?
			end
		end

		@valid
	end

	def invalid?
		!valid?
	end

	def seq
		@seq ||= Bio::Sequence::NA.new contig.seq[(start-1)..(finish-1)]
	end

	def aa_seq(frame)
		seq.translate frame
	end

	def forward?
		direction == '+'
	end

	def direction
		@direction ||= begin
			middle = (sl_mapping.start+sl_mapping.finish)/2.0
			start_dist = (start - middle).abs
			finish_dist = (finish - middle).abs
			start_dist <= finish_dist ? '+' : '-'
		end
	end

	def local_frame
		@local_frame ||= begin
			frames = forward? ? [1, 2, 3] : [4, 5, 6]

			indexes = frames.map { |f| aa_seq(f).index(START_CODON) }
			min = indexes.compact.min

			if min
				min_start_id = indexes.index(min)
				frames[min_start_id]
			end
		end
	end

	def global_frame
		@global_frame ||= begin
			left_idx = forward? ? start : finish

			frame_mapping = {
				'+' => {
							 		0 => { 1 => 3, 2 => 1, 3 => 2 },
							 		1 => { 1 => 1, 2 => 2, 3 => 3 },
							 		2 => { 1 => 2, 2 => 3, 3 => 1 }
							 },
				'-' => {
							 		0 => { 4 => 6, 5 => 4, 6 => 5 },
							 		1 => { 4 => 5, 5 => 6, 6 => 4 },
							 		2 => { 4 => 4, 5 => 5, 6 => 6 }
							 }
			}

			frame_mapping[direction][left_idx%3][local_frame]
		end
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
		@best_blast_hit ||= begin
			blast_hits.sort do |a, b|
				[b.data.data[:evalue].to_f, a.data.data[:pident].to_f] <=> [a.data.data[:evalue].to_f, b.data.data[:pident].to_f]
			end.last
		end
	end

	def sl_mappings
		@sl_mappings ||= begin
			zones = [[start-outer_threshold, start+inner_threshold],
							 [finish-inner_threshold, finish+outer_threshold]]
			contig.sl_mappings.select_intersected(*zones)
		end
	end

	def sl_mapping
		sl_mappings.first
	end

	def to_gff
		color = valid? ? VALID_COLOR : INVALID_COLOR

		if annotated?
			left, right = [gene_start, gene_finish].sort
			f = global_frame
			d = direction
			notes = "ID=#{SecureRandom.hex}"
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

		f = [4,5,6].include?(f) ? (f - 4) : (f - 1)
		[contig.title, :blast, :gene, left, right, '.', d, f, notes].join("\t")
	end

	def normalize!
		@start, @finish = [@start, @finish].sort
	end

	def to_range
		start..finish
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

		hit_cluster.extra_data[:frame] == global_frame
	end

	def inner_threshold
		Settings.annotator.zoi_sl_searching_inner_threshold.to_i
	end

	def outer_threshold
		Settings.annotator.zoi_sl_searching_outer_threshold.to_i
	end
end
