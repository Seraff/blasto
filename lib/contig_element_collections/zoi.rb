module ContigElementCollections
	class Zoi < ContigElementCollection
		attr_accessor :filtered

		def prepare
			@filtered = {}

			normalize
			filter_by_size
			split_policistronic

			log_prepared

			filter_without_blaster
			merge_duplicates
			filter_totally_covered
			filter_intersected

			select(&:invalid?).each do |z|
				@filtered[z.validation_error] ||= []
				@filtered[z.validation_error] << z
			end
			keep_if(&:valid?)

			log_filtered

			# each do |z|
			# 	write_log_to_new_gff 'tmp/annotator/contigs/NODE_1_length_304652_cov_71.8364/a_test.gff',
			# 											 z.best_blast_hit.start,
			# 											 z.best_blast_hit.finish,
			# 											 extra: { 'color' => '#e510ed' }
			# end

			self
		end

		def find_by_id(name)
			detect { |e| e.id.start_with? name }
		end

		def filtered
			self.class.new @filtered.values.flatten
		end

		# protected

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

		def filter_by_size
			small = []

			each do |e|
				small << e if (e.finish - e.start + 1) < Settings.annotator.transcriptome_min_size
			end

			@filtered[:short] = small

			delete_if { |e| small.include? e }
		end

		def split_policistronic
			splitted_zois = []
			zois_for_delition = []

			each do |e|
				if e.polycistronic?
					splitted_zois += e.split_by_polycistronic_cutting_places

					@filtered[:multiple_blaster_groups] ||= []
					@filtered[:multiple_blaster_groups] << e

					zois_for_delition << e
				end
			end

			push(*splitted_zois)
			zois_for_delition.each { |z| delete(z) }
		end

		def filter_without_blaster
			without_blaster = []

			each do |e|
				without_blaster << e if e.hit_clusters.count == 0
			end

			@filtered[:no_hit_clusters] = without_blaster

			delete_if { |e| without_blaster.include? e }
		end

		def merge_duplicates
			elements = {}

			each do |e|
				elements[e.to_range] = e
			end

			vals = elements.values
			keep_if { |e| vals.include? e }
		end

		def filter_totally_covered
			covered = []

			each do |e|
				add_to_covered = false

				other = self.dup
	      other.delete e

	      other.each do |o|
	        if covers? o.to_range, e.to_range
	          add_to_covered = true
	          break
	        end
	      end

	      covered << e if add_to_covered
			end

			@filtered[:totally_covered] = covered

			delete_if { |e| covered.include? e }
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
