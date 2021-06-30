class ReadsStatistics
  attr_reader :reads_file_name

  def initialize(file_name)
    @reads_file_name = file_name
  end

  def average_coverage(contig_name, left, right)
    result = `samtools depth -a -r '#{contig_name}:#{left}-#{right}' #{@reads_file_name} | cut -f3`
    values = result.split("\n").map(&:to_i)
    DescriptiveStatistics::Stats.new(values).median
  end
end
