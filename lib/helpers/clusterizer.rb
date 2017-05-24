require_relative './helpers.rb'

class Clusterizer
	attr_accessor :max_distance, :input_path, :output_path
	attr_reader :sorted, :clustered, :merged

	TMP_DIR = 'tmp/'

	def initialize(input:, output:, max_distance:)
    @input_path = input
    @output_path = output
    @max_distance = max_distance

    @session = SecureRandom.hex
  end

  def perform
    create_tmp_folder
    do_perform
  ensure
    remove_tmp_folder
  end

  def do_perform
  	sort && cluster && merge
  end

  def sort
  	raise NotImplementedError
	end

	def cluster
		raise NotImplementedError
	end

	def merge
		raise NotImplementedError
	end

  protected

  def sorted_path
    change_path @input_path, new_dir: tmp_dir, append: "sorted"
  end

  def clustered_path
    change_path sorted_path, append: "clustered"
  end

  def merged_path
    change_path clustered_path, append: 'merged'
  end

  def tmp_dir
    @tmp_dir ||= Pathname.new(TMP_DIR) + Pathname.new(@session)
  end

  def create_tmp_folder
    dir_path = Pathname.new(tmp_dir)
    dir_path.rmtree if Pathname.new(tmp_dir).exist?
    dir_path.mkpath
  end

  def remove_tmp_folder
    Pathname.new(tmp_dir).rmtree
  end

  def merge_file(input_file, output_file, start_index:, finish_index:, delimiter: "\t", comment_prefix: '#')
  	File.open(output_file, 'w') do |outfile|
	    current_data = nil

	    File.open(input_file, 'r').each do |entry|
	      next if entry.start_with? comment_prefix
	      data = entry.split(delimiter)

	      unless current_data
	        current_data = data
	        next
	      end

	      if current_data[-1] == data[-1]
	        current_data[start_index] = data[start_index].to_i if data[start_index].to_i < current_data[start_index].to_i
	        current_data[finish_index] = data[finish_index].to_i if data[finish_index].to_i > current_data[finish_index].to_i
	      else
	        outfile.puts current_data[0..-2].join(delimiter)
	        current_data = data
	      end
	    end

	    outfile.puts current_data[0..-2].join(delimiter) if current_data
	  end
  end
end

Dir["#{ROOT_PATH}/lib/helpers/clusterizers/*.rb"].each {|file| require file }
