module Clusterizers
  class Bed < Clusterizer
  	START_COL_ID = 1
    FINISH_COL_ID = 2

  	def sort
  		`bedtools sort -i #{@input_path} > #{sorted_path}`
      @sorted = true
  	end

  	def cluster
  		raise 'Not sorted' unless @sorted

  		`bedtools cluster -d #{@max_distance} -i #{sorted_path} > #{clustered_path}`
		  @clustered = true
  	end

  	def merge
  		merge_file clustered_path, output_path, start_index: START_COL_ID, finish_index: FINISH_COL_ID
  	end
  end
end
