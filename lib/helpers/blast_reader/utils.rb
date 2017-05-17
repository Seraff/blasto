class BlastReader
  module Utils
    def back_translate!(target:, extend_borders: false, progress_bar: nil, output_path: nil)
      if output_path
        output_file = File.open(output_path, 'w')
        output_file.puts headers.join(delimiter)
      end

      with_file_rewinded do
        each_hit do |hit|
          hit.extend_borders! target if extend_borders
          hit.back_translate_coords! target

          output_file.puts hit.to_csv if output_file
          progress_bar.increment if progress_bar
        end
      end
    end

    def back_translate_to_gff(output_file, target:, mode:, progress_bar: nil, extend_borders: false)
      raise "Invalid target" unless %w(query subject).include?(target.to_s)
      raise "Invalid mode" unless %w(genome transcriptome).include?(mode.to_s)

      with_file_rewinded do
        each_hit do |hit|
          hit.extend_borders! target if extend_borders
          hit.back_translate_coords! target

          gff = hit.to_gff target, extra_data_keys: [:evalue]

          # additional gff processing
          gff_array = gff.split("\t")

          case mode.to_sym
          when :genome
            gff_array[0].gsub!(/_\d+\z/, '')
          when :transcriptome
            gff_array[0].gsub!(/_length_.+/, '')
          end

          output_file.puts gff_array.join("\t")

          progress_bar.increment if progress_bar
        end
      end
    end

    ## merge hits from one reference contig
    def merge_hits(output_file, target:, max_distance: 128, progress_bar: nil)
      prefix = target.to_s[0] # 'q' or 's'

      sort_by! qseqid: :string, sseqid: :string, "#{prefix}frame" => :digit, "#{prefix}start" => :digit

      merged_hits = []
      merging_group = []

      hits.each do |hit|
        if merging_group.empty?
          merging_group << hit

        elsif hit.needs_to_merge_with? merging_group.last, target: target, max_distance: max_distance
          merging_group << hit

        else
          # merge all hits in group
          if merging_group.count == 1
            merged_hits << merging_group.first
          else
            merged_hits << merge_hit_group(merging_group, target: target)
          end
          merging_group = [hit]
        end
        progress_bar.increment if progress_bar
      end

      output_file.puts headers.join(delimiter)

      merged_hits.each do |hit|
        output_file.puts hit.to_csv(delimiter: delimiter)
      end

      output_file.close
    end

    def sort_by!(keys, ouput_path: nil)
      keys = [keys] unless keys.is_a? Array
      indexes = keys.map { |k, type| [@headers.index(k.to_s), type] }.to_h
      raise 'Incorrect keys' if indexes.any?(&:nil?)

      indexes = indexes.map { |i, type| [i+1, type] }.to_h
      index_str = indexes.map { |i, type| "-k #{i},#{i}#{ type == :digit ? 'n' : '' }" }.join ' '
      new_path = ouput_path || change_path(file.path, append: 'tmp')

      `(head -n 1 #{file.path} && (tail -n +2 #{file.path} | sort -t#{delimiter} #{index_str})) > #{new_path}`
      `mv #{new_path} #{file.path}` unless ouput_path

      reopen ouput_path
    end
  end
end
