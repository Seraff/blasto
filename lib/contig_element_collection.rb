require_relative 'contig_element'

class ContigElementCollection < Array
	attr_accessor :collection

	def select_intersected(*coords)
		min_start = coords.map{ |e| e[0] }.min
		max_finish = coords.map{ |e| e[1] }.max

		coll = select { |e| e.start <= max_finish }.select { |e| e.finish >= min_start }

		result = coll.select do |e|
			coords.select { |c| ![c[0]..c[1], (e.start)..(e.finish)].intersection.nil? }.any?
		end

		self.class.new result
	end

	def contig=(contig)
		@contig = contig
	end

	# smells very bad
	def covers?(a, b)
    [a, b].intersection == (b.first..b.last)
  end
end

Dir["#{ROOT_PATH}/lib/contig_element_collections/*.rb"].each {|file| require file }
