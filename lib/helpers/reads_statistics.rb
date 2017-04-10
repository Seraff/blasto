class ReadsStatistics
  attr_reader :reads_file_name

  def initialize(file_name)
    @reads_file_name = file_name
  end

  def median_coverage(contig_name, left, right)
    result = `samtools depth -r '#{contig_name}:#{left}-#{right}' #{@reads_file_name}`
    values = result.split("\n").map { |r| r.split(/\s+/)[-1].to_i }

    DescriptiveStatistics::Stats.new(values).median
  end
end
