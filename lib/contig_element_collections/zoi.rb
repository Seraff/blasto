module ContigElementCollections
	class Zoi < ContigElementCollection
		attr_accessor :filtered

		def prepare
			normalize
			split_polycistronic
			log_prepared
			merge_duplicates

			self
		end

		def find_by_id(name)
			detect { |e| e.id.start_with? name }
		end

		protected

		def normalize
			each(&:normalize!)
		end

		#TODO: what to do? we do not use it
		def merge_by_best_blast_hit
			to_remove = []
			merged = []

			group_by(&:best_blast_hit).each do |hit, group|
				next if hit.nil?

				if group.count > 1
					group.each { |e| e.make_invalid! reason: :merged }

					start = group.sort_by { |e| e.start }.first.start
					finish = group.sort_by { |e| e.finish }.last.finish
					merged << ContigElements::Zoi.new(group.first.contig,
						                                start,
					                                  finish,
					                                  group.first.raw_gff)
				end
			end

			push(*merged)
		end

		def split_polycistronic
			splitted_zois = []

			each do |e|
				if e.polycistronic?
					splitted_zois += e.split_by_polycistronic_cutting_places
					e.make_invalid! reason: :splitted
				end
			end

			push(*splitted_zois)
		end

		def merge_duplicates
			elements = {}

			each do |e|
				elements[e.to_range] = e
			end

			vals = elements.values
			keep_if { |e| vals.include? e }
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
