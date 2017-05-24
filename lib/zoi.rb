class Zoi
	DIRECTIONS = ['+', '-'].freeze
	FRAMES = (1..6).to_a.freeze

	attr_reader :contig, :start, :finish, :direction, :frame

	def initialize(contig, start, finish)
		@valid = true

		@contig = contig
		@start = start
		@finish = finish
	end

	def valid?
		hit_clusters.one? && !sl_mapping.nil?
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

	def frame
		@frame ||= begin
			frames = forward? ? [1, 2, 3] : [4, 5, 6]

			indexes = frames.map { |f| aa_seq(f).index('M') }
			min_start_id = indexes.index(indexes.compact.min)
			local_frame = frames[min_start_id]

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
		if forward?
			# TODO
		else
			finish
		end
	end

	def gene_finish
		if forward?
			finish
		else
			# TODO
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
			contig.sl_mappings.select_intersected([start-threshold, start+threshold], [finish-threshold, finish+threshold])
		end
	end

	def sl_mapping
		sl_mappings.first
	end

	def to_gff
		f = [4,5,6].include?(frame) ? (frame - 4) : (frame - 1)
		[contig.title, :blast, :gene, start, finish, '.', direction, f, "ID=#{SecureRandom.hex}_#{frame}"].join("\t")
	end

	protected

	def threshold
		Settings.annotator.zoi_sl_searching_threshold.to_i
	end
end
