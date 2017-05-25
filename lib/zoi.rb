class Zoi
	DIRECTIONS = ['+', '-'].freeze
	FRAMES = (1..6).to_a.freeze
	START_CODON = 'M'

	attr_reader :contig, :start, :finish, :direction, :frame, :invalidity_reason, :extra_data

	def initialize(contig, start, finish, extra_data: {})
		@contig = contig
		@start = start
		@finish = finish
		@invalidity_reason = nil
		@extra_data = extra_data
	end

	def valid?
		@invalidity_reason = nil

		if !hit_clusters.one?
			@invalidity_reason = :hit_clusters_not_one
		elsif sl_mapping.nil?
			@invalidity_reason = :has_no_sl_mappings
		elsif local_frame.nil?
			@invalidity_reason = :cannot_detect_frame
		elsif !in_bh_cluster_frame?
			@invalidity_reason = :hit_cluster_has_another_frame
		end

		if @invalidity_reason
			BadTranscriptsLogger.add_to_bin contig.title, raw_gff: extra_data[:raw_gff],
																										reason: @invalidity_reason
		end

		@invalidity_reason.nil?
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

	def gene_start
		idx = aa_seq(local_frame).index(START_CODON)

		if forward?
			start+idx*3+local_frame-1
		else
			finish-idx*3-local_frame+4
		end
	end

	def gene_finish
		if forward?
			finish
		else
			start
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
		left, right = [gene_start, gene_finish].sort
		f = [4,5,6].include?(global_frame) ? (global_frame - 4) : (global_frame - 1)
		[contig.title, :blast, :gene, left, right, '.', direction, f, "ID=#{SecureRandom.hex}_#{global_frame};color=#009900"].join("\t")
	end

	protected

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
