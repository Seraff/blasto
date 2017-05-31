module ContigElements
	class Basic < ContigElement
		attr_accessor :data, :extra_data

		def initialize(contig, start, finish, data, extra_data: {})
			@contig = contig
			@start = [start, finish].map(&:to_i).min
			@finish = [start, finish].map(&:to_i).max
			@data = data
			@extra_data = extra_data
		end
	end
end
