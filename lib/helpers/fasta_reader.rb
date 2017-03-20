require 'bio'

class FastaReader < Bio::FlatFile
  def initialize(file_name)
    super Bio::FastaFormat, File.open(file_name)
  end

  def gc_content
    rewind

    gc_count = 0
    at_count = 0

    each do |e|
      gc_count += e.seq.split('').count { |e| %w(G C).include? e }
      at_count += e.seq.split('').count { |e| %w(A T U).include? e }
    end

    rewind

    (gc_count.to_f/(gc_count.to_f+at_count.to_f))
  end

  # smells bad: what about reverse strand?..
  def codon_usage
    rewind

    total_usage = {}

    each do |e|
      Bio::Sequence::NA.new(e.seq).codon_usage.each do |c, count|
        aa = Bio::Sequence::NA.new(c).translate
        total_usage[aa] ||= {}
        total_usage[aa][c.upcase] ||= 0
        total_usage[aa][c.upcase] += count
      end
    end

    rewind

    total_usage
  end

  def to_h
    @to_h ||= begin
      rewind
      result = map { |c| [c.entry_id, c.seq] }.to_h
      rewind
      result
    end
  end
end
