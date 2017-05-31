module ContigElements
	class Zoi < ContigElement
		module Polycistronic
			def polycistronic?
			end

			def multiple_blaster_groups?
				polycistronic_sl_mappings.any?
			end

			def split_by_polycistronic_sl_mappings
				return false unless multiple_blaster_groups?
				# TODO
			end

			def polycistronic_sl_mappings
				@polycistronic_sl_mappings ||= begin
					return [] if hit_clusters.count <= 1 || all_sl_mappings.empty?

					groups = group_blasters_by_intersections
					return [] if groups.count <= 1

					# ---[***]------[***]---
					# ----------SL----------

					sls_for_cutting = []

					groups.each_with_index do |group, i|
						next_group = groups[i+1]
						break if next_group.nil?

						left = group.max { |e| e.finish }.finish - Settings.annotator.polycistronic_sl_threshold
						right = next_group.min { |e| e.start }.start + Settings.annotator.polycistronic_sl_threshold

						sls = all_sl_mappings.dup.select_intersected([left, right])

						sls_for_cutting << sls.first if sls.any?
					end

					sls_for_cutting
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
		end
	end
end
