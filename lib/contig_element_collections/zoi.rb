module ContigElementCollections
	class Zoi < ContigElementCollection
		attr_accessor :filtered

		def prepare
			normalize
			check_validation_defection
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

		def check_validation_defection
			each(&:validate)
			each(&:check_defection)
		end

		def split_polycistronic
			splitted_zois = []

			each do |e|
				next unless e.valid?

				if e.polycistronic?
					splitted_zois += e.polycistronic_subzois
					e.make_invalid! reason: :splitted
				end
			end

			splitted_zois.each do |z|
				z.make_defective! reason: :produced_by_splitting
				z.make_defective! reason: :fused_genes if z.blast_hits.count > 1
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
