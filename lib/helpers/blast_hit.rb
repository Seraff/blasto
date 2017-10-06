require_relative './blast_hit_data.rb'
require_relative './blast_hit/merging.rb'
require_relative './blast_hit/accessors.rb'
require_relative './blast_hit/gff.rb'

class BlastHit
  include Merging
  include Accessors
  include Gff

  TARGET_KEYS = {
    query:   { id: :qseqid, start: :qstart, finish: :qend, opposite_id: :sseqid, frame: :qframe, len: :qlen },
    subject: { id: :sseqid, start: :sstart, finish: :send, opposite_id: :qseqid, frame: :sframe, len: :slen }
  }

  attr_accessor :data
  attr_reader :headers, :extended

  def initialize(headers, data)
    @headers = headers
    @data = BlastHitData.new data
    @extended = false
  end

  def merging_gaps
    @merging_gaps ||= []
    @merging_gaps
  end

  def merged?
    merging_gaps.any?
  end

  def unextended_data
    @unextended_data ||= {}
    @unextended_data
  end

  def extended?
    unextended_data.any?
  end

  # everything should be in AA coords!
  def extend_borders!(target)
    keys = detect_keys target
    opposite_keys = opposite_keys target

    left, right = [data[keys[:start]], data[keys[:finish]]].sort
    op_left, op_right = [data[opposite_keys[:start]], data[opposite_keys[:finish]]].sort
    op_len = data[opposite_keys[:len]]

    left_shift = op_left
    right_shift = op_len-op_right

    new_left = left - left_shift
    new_left = 1 if new_left < 1

    new_right = right + right_shift
    new_right = data[keys[:len]] if new_right > data[keys[:len]]

    @extended = true
    @unextended_data ||= {}
    @unextended_data[keys[:start]] = data[keys[:start]]
    @unextended_data[keys[:finish]] = data[keys[:finish]]

    data[keys[:start]] = data[keys[:start]] < data[keys[:finish]] ? new_left : new_right
    data[keys[:finish]] = data[keys[:start]] < data[keys[:finish]] ? new_right : new_left
  end

  def back_translate_coords!(target)
    keys = detect_keys target
    start = data[keys[:start]]
    finish = data[keys[:finish]]

    data[keys[:start]], data[keys[:finish]] = back_translate_coord_pair(target, start, finish)
    data[keys[:len]] = detect_nalen target

    @unextended_data ||= {}
    if @unextended_data.any?
      unext_start = @unextended_data[keys[:start]]
      unext_finish = @unextended_data[keys[:finish]]
      unext_start, unext_finish = back_translate_coord_pair(target, unext_start, unext_finish)
      @unextended_data[keys[:start]] = unext_start
      @unextended_data[keys[:finish]] = unext_finish
    end

    @merging_gaps ||= []
    new_gaps = []
    @merging_gaps.each do |gap|
      new_gaps << back_translate_coord_pair(target, gap[0], gap[1])
    end
    @merging_gaps = new_gaps
  end

  def detect_keys(target)
    raise "Wrong target" if !target || !TARGET_KEYS.keys.include?(target.to_sym)
    target = target.to_sym
    TARGET_KEYS[target]
  end

  def attr_by_target(attr_name, target)
    data[TARGET_KEYS[target.to_sym][attr_name.to_sym]]
  end

  def assign_attr_by_target(attr_name, value, target)
    @data[TARGET_KEYS[target.to_sym][attr_name.to_sym]] = value
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

  def organism(target)
    org = data[detect_keys(target)[:id]]
    return nil unless org.include?('|')
    org.split('|').first
  end

  def evalue
    @evalue ||= BigDecimal.new(data[:evalue].to_s)
  end

  protected

  def back_translate_coord_pair(target, a, b)
    frame = detect_frame target
    nalen = detect_nalen target

    transposed_frame = [4,5,6].include?(frame) ? frame - 3 : frame

    new_a = 3*a+transposed_frame-3
    new_b = 3*b+transposed_frame-1

    if [4,5,6].include? frame
      new_a, new_b = detect_reversed_coords new_a, new_b, nalen
    end

    [new_a, new_b]
  end

  def detect_reversed_coords(start, finish, len)
    [start, finish].map { |x| len+1-x }
  end
end
