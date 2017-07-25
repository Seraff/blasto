module Clusterizers
  class Gff < Clusterizer
    GFF_HEADER = '##gff-version 3'

    START_COL_ID = 3
    FINISH_COL_ID = 4

    def do_perform
      sort && cluster && merge && join_files
    end

    def sort
      `bedtools sort -i #{@input_path} > #{sorted_path}`
      @sorted = true
    end

    def cluster
      raise 'Not sorted' unless @sorted

      FRAMES.each do |frame|
        command = ["cat #{sorted_path}"]
        command << "grep -P \"#{regexp_for_frame(frame)}\""
        command << "bedtools cluster -d #{@max_distance}"
        command = command.join(' | ')
        `#{command} > #{clustered_path(frame)}`
      end

      @clustered = true
    end

    def merge
      raise 'Not sorted or not clustered' if !@sorted || !@clustered

      FRAMES.each do |frame|
        merge_file clustered_path(frame), merged_path(frame), start_index: START_COL_ID, finish_index: FINISH_COL_ID
      end
    end

    def join_files
      paths = FRAMES.map { |f| merged_path(f) }
      raise 'Merged tmp files does not present' if paths.any? { |p| not Pathname.new(p).exist? }

      `printf "#{GFF_HEADER}\n" > #{output_path}`
      `cat #{paths.join(' ')} >> #{output_path}`
    end

    protected

    def regexp_for_frame(frame)
      ['(.+?\t){3}[+-]\t(', frame, ')'].join
    end

    def clustered_path(frame)
      change_path sorted_path, append: "clustered_#{frame}"
    end

    def merged_path(frame)
      change_path clustered_path(frame), append: 'merged'
    end
  end
end
