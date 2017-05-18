require_relative './helpers.rb'

class GffClusterizer
  attr_accessor :max_distance, :input_path, :output_path
  attr_reader :sorted, :clustered, :merged

  TMP_DIR = 'tmp/'
  FRAMES = (1..6).to_a
  GFF_HEADER = '##gff-version 3'

  def initialize(input:, output:, max_distance:)
    @input_path = input
    @output_path = output
    @max_distance = max_distance

    @session = SecureRandom.hex
  end

  # the main action
  def cluster_and_merge
    create_tmp_folder
    sort && cluster && merge && join_files
  ensure
    remove_tmp_folder
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
      outfile = File.open(merged_path(frame), 'w')

      current_data = nil

      File.open(clustered_path(frame), 'r').each do |entry|
        next if entry.start_with? '#'
        data = entry.split("\t")

        unless current_data
          current_data = data
          next
        end

        if current_data[-1] == data[-1]
          current_data[3] = data[3].to_i if data[3].to_i < current_data[3].to_i
          current_data[4] = data[4].to_i if data[4].to_i > current_data[4].to_i
        else
          outfile.puts current_data[0..-2].join("\t")
          current_data = data
        end
      end

      outfile.puts current_data[0..-2].join("\t") if current_data

      outfile.close
    end
  end

  def join_files
    paths = FRAMES.map { |f| merged_path(f) }
    raise 'Merged tmp files does not present' if paths.any? { |p| not Pathname.new(p).exist? }

    `printf "#{GFF_HEADER}\n" > #{output_path}`
    `cat #{paths.join(' ')} >> #{output_path}`
  end

  protected

  def sorted_path
    change_path @input_path, new_dir: tmp_dir, append: "sorted"
  end

  def clustered_path(frame)
    change_path sorted_path, append: "clustered_#{frame}"
  end

  def merged_path(frame)
    change_path clustered_path(frame), append: 'merged'
  end

  def tmp_dir
    @tmp_dir ||= Pathname.new(TMP_DIR) + Pathname.new(SecureRandom.hex)
  end

  def create_tmp_folder
    dir_path = Pathname.new(tmp_dir)
    dir_path.rmtree if Pathname.new(tmp_dir).exist?
    dir_path.mkpath
  end

  def remove_tmp_folder
    Pathname.new(tmp_dir).rmtree
  end

  def regexp_for_frame(frame)
    ['(.+?\t){3}[+-]\t(', frame, ')'].join
  end
end
