require_relative 'contig_element'
require_relative 'zoi'

class ContigElementCollection < Array
	attr_accessor :collection

	def select_intersected(*coords)
		min_start = coords.map{ |e| e[0] }.min
		max_finish = coords.map{ |e| e[1] }.max

		coll = select { |e| e.start <= max_finish }.select { |e| e.finish >= min_start }

		coll.select do |e|
			coords.select { |c| ![c[0]..c[1], (e.start)..(e.finish)].intersection.nil? }.any?
		end
	end
end
