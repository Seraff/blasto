require_relative './blast_hit_data.rb'

class BlastHit
  TARGET_KEYS = {
    query:   { id: :qseqid, start: :qstart, finish: :qend, opposite_id: :sseqid },
    subject: { id: :sseqid, start: :sstart, finish: :send, opposite_id: :qseqid }
  }

  attr_accessor :data

  def initialize(data)
    @data = BlastHitData.new data
  end

  def to_gff(target)
    keys = detect_keys target

    required_fields = [keys[:id], keys[:start], keys[:finish]]
    required_fields.each do |f|
      raise 'Cannot convert: not enough fields' unless data.key?(f)
    end

    start = data[keys[:start]]
    finish = data[keys[:finish]]

    strand = '+'
    if start > finish
      strand = '-'
      tmp = start
      start = finish
      finish = tmp
    end

    id = data[keys[:id]]
    start = 1 if start.zero?
    frame = detect_frame target

    note = "ID=#{data[keys[:opposite_id]]}_#{(1..16).to_a.map{ (0..9).to_a.sample }.join}_#{frame}"

    [id, 'blast', 'gene', start, finish, '.', strand, frame, note].join("\t")
  end

  def back_translate_coords!(target)
    keys = detect_keys target

    frame = detect_frame target
    nlen = detect_nalen target
    start = data[keys[:start]]
    finish = data[keys[:finish]]

    transposed_frame = [4,5,6].include?(frame) ? frame - 3 : frame
    new_start = 3*start+transposed_frame-3
    new_finish = 3*finish+transposed_frame-1

    if [4,5,6].include? frame
      new_start, new_finish = detect_reversed_coords new_start, new_finish, nlen
    end

    data[keys[:start]] = new_start
    data[keys[:finish]] = new_finish
  end

  def detect_keys(target)
    target = target.to_sym
    raise "Wrong target" unless TARGET_KEYS.keys.include? target
    TARGET_KEYS[target]
  end

  def detect_frame(target)
    data[detect_keys(target)[:id]].split('_')[-1].to_i
  end

  def detect_nalen(target)
    data[detect_keys(target)[:id]].match(/length_(?<len>\d+)/)[:len].to_i
  end

  protected

  def detect_reversed_coords(start, finish, len)
    [start, finish].map { |x| len+1-x }
  end
end
