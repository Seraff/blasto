require_relative './blast_hit_data.rb'

class BlastHit
  attr_accessor :data

  def initialize(data)
    @data = BlastHitData.new data
  end

  def to_gff(target: :query, frame: '0')
    case target
    when :query
      id_key = :qseqid
      start_key = :qstart
      finish_key = :qend
      note_key = :sseqid
    when :subject
      id_key = :sseqid
      start_key = :sstart
      finish_key = :send
      note_key = :qseqid
    end
    required_fields = [id_key, start_key, finish_key]

    required_fields.each do |f|
      raise 'Cannot convert: not enough fields' unless data.key?(f)
    end

    start = data[start_key]
    finish = data[finish_key]

    strand = '+'
    if start > finish
      strand = '-'
      tmp = start
      start = finish
      finish = tmp
    end

    id = data[id_key]
    start = 1 if start.zero?

    note = "ID=#{data[note_key]}_#{(1..16).to_a.map{ (0..9).to_a.sample }.join}_#{frame}"

    [id, 'blast', 'gene', start, finish, '.', strand, frame, note].join("\t")
  end
end
