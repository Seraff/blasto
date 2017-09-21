class BlastReader
  module Utils
    def back_translate!(target:, progress_bar: nil, output_path: nil)
      perform_with_each_hit(progress_bar: progress_bar, output_path: output_path) do |hit|
        hit.back_translate_coords! target
      end
    end

    def extend_borders!(target:, progress_bar: nil, output_path: nil)
      perform_with_each_hit(progress_bar: progress_bar, output_path: output_path) do |hit|
        hit.extend_borders! target
      end
    end

    def merge_hits!(target:, output_path: nil, max_distance: 256, progress_bar: nil)
      prefix = target.to_s[0] # 'q' or 's'
      sort_by! qseqid: :string, sseqid: :string, "#{prefix}frame" => :digit, "#{prefix}start" => :digit

      cache_hits unless hits_cached?

      merged_hits = []
      merging_group = []

      each_hit do |hit|
        if merging_group.empty?
          merging_group << hit

        elsif hit.needs_to_merge_with? merging_group.last, target: target, max_distance: max_distance
          merging_group << hit

        else
          # merge all hits in a group
          if merging_group.count == 1
            merged_hits << merging_group.first
          else
            merged_hits << merge_hit_group(merging_group, target: target)
          end
          merging_group = [hit]
        end

        progress_bar.increment if progress_bar
      end

      merged_hits << (merging_group.one? ? merging_group.first : merge_hit_group(merging_group, target: target))

      self.hits = merged_hits

      if output_path
        output_file = File.open(output_path, 'w')
        output_file.puts headers.join(delimiter)

        merged_hits.each do |hit|
          output_file.puts hit.to_csv(delimiter: delimiter)
        end

        output_file.close
      end
    end

    def merge_hit_group(hits, target:)
      # puts "Merging hits group: query #{hits.first.data[:qseqid]}, subject #{hits.first.data[:sseqid]}"
      hit = hits.first
      hits = hits[1..-1]

      hits.each do |h|
        hit.merge_with h, target: target
      end

      hit
    end

    def sort_by!(keys, ouput_path: nil)
      if hits_cached?
        @hits.sort_by! { |h| keys.keys.map { |k| h.data[k] } }
      else
        indexes = keys.map { |k, type| [@headers.index(k.to_s), type] }.to_h

        raise 'Incorrect keys' if indexes.keys.any?(&:nil?)

        indexes = indexes.map { |i, type| [i+1, type] }.to_h
        index_str = indexes.map { |i, type| "-k #{i},#{i}#{ type == :digit ? 'n' : '' }" }.join ' '
        new_path = ouput_path || change_path(file.path, append: 'tmp')

        `(head -n 1 #{file.path} && (tail -n +2 #{file.path} | sort -t#{delimiter} #{index_str})) > #{new_path}`
        `mv #{new_path} #{file.path}` unless ouput_path

        reopen ouput_path
      end
    end

    protected

    def perform_with_each_hit(progress_bar: nil, output_path: nil)
      if output_path
        output_file = File.open(output_path, 'w')
        output_file.puts headers.join(delimiter)
      end

      with_file_rewinded do
        each_hit do |hit|
          yield hit

          output_file.puts hit.to_csv if output_file
          progress_bar.increment if progress_bar
        end
      end

      output_file.close if output_path
    end
  end
end
