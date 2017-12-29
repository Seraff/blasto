module ContigElementCollections
	class Zoi < ContigElementCollection
		attr_accessor :filtered

		def prepare
			@filtered = {}

			normalize
			split_polycistronic

			log_prepared

			merge_duplicates
			filter_totally_covered
			filter_intersected

			select(&:invalid?).each do |z|
				@filtered[z.validation_error] ||= []
				@filtered[z.validation_error] << z
			end
			keep_if(&:valid?)

			log_filtered

			self
		end

		def find_by_id(name)
			detect { |e| e.id.start_with? name }
		end

		def filtered
			self.class.new @filtered.values.flatten
		end

		protected

		def normalize
			each(&:normalize!)
		end

		def merge_by_best_blast_hit
			to_remove = []
			merged = []

			group_by(&:best_blast_hit).each do |hit, group|
				next if hit.nil?

				if group.count > 1
					to_remove += group

					start = group.sort_by { |e| e.start }.first.start
					finish = group.sort_by { |e| e.finish }.last.finish
					merged << ContigElements::Zoi.new(group.first.contig,
						                                start,
					                                  finish,
					                                  group.first.raw_gff)
				end
			end

			@filtered[:merged] = to_remove
			delete_if { |e| to_remove.include? e }
			push(*merged)
		end

		def split_polycistronic
			splitted_zois = []
			zois_for_delition = []

			each do |e|
				if e.polycistronic?
					splitted_zois += e.split_by_polycistronic_cutting_places

					zois_for_delition << e
				end
			end

			push(*splitted_zois)
			zois_for_delition.each { |z| delete(z) }
		end

		def merge_duplicates
			elements = {}

			each do |e|
				elements[e.to_range] = e
			end

			vals = elements.values
			keep_if { |e| vals.include? e }
		end

		def filter_intersected
			intersected = []

			each do |me|
	      next if intersected.include? me

	      without_me = self.dup
	      without_me.delete(me)

	      intersection_found = false

	      without_me.each do |other|
	        if intersects? other.to_range, me.to_range
	          intersection_found = true
	          intersected += [me, other]
	          break
	        end
	      end
	    end

	    @filtered[:intersected] = intersected

	    delete_if { |e| intersected.include? e }
		end

		def intersects?(a, b)
	  	![a, b].intersection.nil?
		end

	  def covers?(a, b)
	    [a, b].intersection == (b.first..b.last)
	  end

	  def log_filtered
	  	@filtered.each do |reason, group|
	  		group.each do |e|
	  			e.make_invalid! reason: reason
	  			BadTranscriptsLogger.add_to_bin e
	  		end
	  	end
	  end

	  def log_prepared
			return unless Preparer.contig_folder_path(@contig.title).exist?

			File.open(Preparer.contig_folder_path(@contig.title, filename: 'prepared_zois.gff'), 'w') do |f|
				each do |z|
					f.puts z.to_gff_as_is
				end
			end
	  end
	end
end
