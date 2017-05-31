class ContigElement
	attr_accessor :start, :finish, :contig

	def covers?(other)
		IntervalsHelper.covers? (start..finish), (other.start..other.finish)
	end

	def intersects?(other)
		IntervalsHelper.intersects? start..finish, other.start..other.finish
	end

	def seq
		@seq ||= Bio::Sequence::NA.new contig.seq[(start-1)..(finish-1)]
	end

	def aa_seq(frame)
		seq.translate frame
	end

	def na_len
		seq.length
	end
end

Dir["#{ROOT_PATH}/lib/contig_elements/*.rb"].each {|file| require file }
