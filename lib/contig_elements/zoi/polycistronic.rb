module ContigElements
	class Zoi < ContigElement
		module Polycistronic
			class CuttingPlace
				attr_accessor :coord, :sl

				def initialize(coord, sl=nil)
					@coord = coord
					@sl = sl
				end

				def has_sl?
					!@sl.nil?
				end
			end

			def polycistronic?
				polycistronic_cutting_places.any?
			end

			def split_by_polycistronic_cutting_places
				return false unless polycistronic?

				new_zois = []
				last_coord = start

				# puts
				# puts "#{contig.title}"
				# puts "#{polycistronic_cutting_places.map{|a| [a.start, a.finish].join('-')}}, #{start} - #{finish}"

				polycistronic_cutting_places.each do |place|
					new_start = last_coord == start ? start : last_coord + 1
					new_zois << get_subzoi(new_start, place.coord)

					last_coord = place.coord
				end

				new_zois << get_subzoi(last_coord+1, finish)

				new_zois
			end

			def polycistronic_cutting_places
				@polycistronic_cutting_places ||= begin
					return [] if hit_clusters.count <= 1

					groups = group_blasters_by_intersections
					return [] if groups.count <= 1

					# ---[***]------[***]---
					# ----------SL----------

					cutting_places = []

					groups.each_with_index do |group, i|
						next_group = groups[i+1]
						break if next_group.nil?

						left = group.max { |e| e.finish }.finish - Settings.annotator.polycistronic_sl_threshold
						right = next_group.min { |e| e.start }.start + Settings.annotator.polycistronic_sl_threshold

						sl = all_sl_mappings.dup.select_intersected([left, right]).first

						if sl
							cutting_places << CuttingPlace.new((sl.start+sl.finish)/2, sl)
						else
							cutting_places << CuttingPlace.new((left+right)/2)
						end
					end

					cutting_places
				end
			end

			def group_blasters_by_intersections
			  groups = []
			  sorted = hit_clusters.dup.sort_by { |h| h.start }

			  current_group = []
			  sorted.each do |e|
			    found = false

			    if current_group.empty?
			      current_group << e
			      next
			    end

			    current_group.each do |ge|
			      if e.intersects?(ge)
			        current_group << e
			        found = true
			        break
			      end
			    end

			    unless found
			      groups << current_group
			      current_group = [e]
			    end
			  end

			  groups << current_group

			  groups
			end

			def get_subzoi(start_coord, finish_coord)
				raise "Wrong coords (#{start_coord}, #{finish_coord})" if start_coord > finish_coord

				start_coord = [start, start_coord].sort.last
				finish_coord = [finish, finish_coord].sort.first

				obj = self.class.new(contig, start_coord, finish_coord, raw_gff)
				obj.make_defective! reason: :splitted
				obj
			end
		end
	end
end
