class ContigElement
	attr_accessor :start, :finish, :data, :extra_data

	def initialize(start, finish, data, extra_data: {})
		@start = [start, finish].map(&:to_i).min
		@finish = [start, finish].map(&:to_i).max
		@data = data
		@extra_data = extra_data
	end

	def totally_covers?(other)
		[(start..finish), (other.start..other.finish)].intersection == (other.start..other.finish)
	end
end
