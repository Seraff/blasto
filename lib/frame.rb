class Frame
  attr_reader :seq, :start, :finish

  def initialize(seq, start, finish)
    raise 'Wrong seq format' if !seq.is_a?(Bio::Sequence::AA) && !seq.is_a?(Bio::Sequence::NA)

    @seq = seq
    @start = start
    @finish = finish
  end
end
