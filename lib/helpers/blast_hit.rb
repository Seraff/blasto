require_relative './blast_hit_data.rb'
require_relative './blast_hit/merging.rb'

class BlastHit
  include Merging

  TARGET_KEYS = {
    query:   { id: :qseqid, start: :qstart, finish: :qend, opposite_id: :sseqid, frame: :qframe, len: :qlen },
    subject: { id: :sseqid, start: :sstart, finish: :send, opposite_id: :qseqid, frame: :sframe, len: :slen }
  }

  attr_accessor :data

  def initialize(headers, data)
    @headers = headers
    @data = BlastHitData.new data
  end

  def to_gff(target, extra_data_keys: [])
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

    if extra_data_keys.any?
      extra_data = extra_data_keys.select { |k| data[k] }.map { |k| "#{k}: #{data[k]}" }.join(', ')
      note += ";Note=#{extra_data}"
    end

    [id, 'blast', 'gene', start, finish, '.', strand, frame, note].join("\t")
  end

  # everything should be in AA coords!
  def extend_borders!(target)
    keys = detect_keys target
    opposite_keys = opposite_keys target

    left_shift = data[opposite_keys[:start]]
    right_shift = data[opposite_keys[:len]]-data[opposite_keys[:finish]]

    new_start = data[keys[:start]] - left_shift
    new_start = 1 if new_start < 1

    new_finish = data[keys[:finish]] + right_shift
    new_finish = data[keys[:len]] if new_finish > data[keys[:len]]

    data[keys[:start]] = new_start
    data[keys[:finish]] = new_finish
  end

  def back_translate_coords!(target)
    keys = detect_keys target

    frame = detect_frame target
    nalen = detect_nalen target
    start = data[keys[:start]]
    finish = data[keys[:finish]]

    transposed_frame = [4,5,6].include?(frame) ? frame - 3 : frame
    new_start = 3*start+transposed_frame-3
    new_finish = 3*finish+transposed_frame-1

    if [4,5,6].include? frame
      new_start, new_finish = detect_reversed_coords new_start, new_finish, nalen
    end

    data[keys[:start]] = new_start
    data[keys[:finish]] = new_finish
    data[keys[:len]] = nalen
  end

  def detect_keys(target)
    raise "Wrong target" if !target || !TARGET_KEYS.keys.include?(target.to_sym)
    target = target.to_sym
    TARGET_KEYS[target]
  end

  def opposite_keys(target)
    detect_keys opposite_target(target)
  end

  def opposite_target(target)
    case target.to_sym
    when :query
      :subject
    when :subject
      :query
    end
  end

  def detect_frame(target)
    data[detect_keys(target)[:frame]] || data[detect_keys(target)[:id]].split('_')[-1].to_i
  end

  def detect_nalen(target)
    data[detect_keys(target)[:id]].match(/length_(?<len>\d+)/)[:len].to_i
  end

  def to_csv(delimiter: ',')
    @headers.map { |h| "#{data[h]}" }.join delimiter
  end

  protected

  def detect_reversed_coords(start, finish, len)
    [start, finish].map { |x| len+1-x }
  end
end
