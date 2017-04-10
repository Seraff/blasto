class FrameCollection
  attr_reader :frames

  def initialize
    @frames = []
  end

  def push(seq, start, finish)
    @frames << Frame.new(seq, start, finish)
  end

  def each
    @frames.each do |f|
      yield f
    end
  end
end

