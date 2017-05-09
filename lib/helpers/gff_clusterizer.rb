require_relative './helpers.rb'

class GffClusterizer
  attr_accessor :max_distance, :input_path, :output_path
  attr_reader :sorted, :clustered, :merged

  TMP_DIR = 'tmp/'

  def initialize(input:, output:, max_distance:)
    @input_path = input
    @output_path = output
    @max_distance = max_distance

    @session = SecureRandom.hex
  end

  def sort
    `bedtools sort -i #{@input_path} > #{sorted_path}`
    @sorted = true
  end

  def cluster
    raise 'Not sorted' unless @sorted
    `bedtools cluster -s -d #{@max_distance} -i #{sorted_path} > #{clustered_path}`
    @clustered = true
  end

  def merge
    raise 'Not sorted or not clustered' if !@sorted || !@clustered

    outfile = File.open(output_path, 'w')
    outfile.puts "##gff-version 3"

    current_data = nil

    File.open(clustered_path, 'r').each do |entry|
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

    outfile.puts current_data[0..-2].join("\t")

    outfile.close
  end

  # the main action
  def cluster_and_merge
    sort && cluster && merge
    remove_tmp_files
  end

  protected

  def sorted_path
    change_path(@input_path, new_dir: TMP_DIR, append: "#{@session}_sorted")
  end

  def clustered_path
    change_path(sorted_path, new_dir: TMP_DIR, append: "clustered")
  end

  def remove_tmp_files
    `rm #{sorted_path} #{clustered_path}`
  end
end
