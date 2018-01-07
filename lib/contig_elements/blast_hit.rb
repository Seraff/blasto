module ContigElements
  class BlastHit < Basic
    def frame
      @frame ||= data.frame(contig.target).to_i
    end

    def organism
      @organism ||= data.data[:qseqid].split(/_|\|/).first
    end

    def gene
      @gene ||= data.data[:qseqid].split(/_|\|/)[1..-1].join
    end

    def direction
      if [1, 2, 3].include?(frame)
        '+'
      elsif [4, 5, 6].include?(frame)
        '-'
      else
        raise "Wrong frame for #{self.class.inspect}##{__method__}: #{frame}"
      end
    end

    def forward?
      direction == '+'
    end

    def reverse?
      !forward
    end

    def extended_start
      cache_extended_borders
      @extended_start
    end

    def extended_finish
      cache_extended_borders
      @extended_finish
    end

    def begin
      forward? ? start : finish
    end

    def end
      forward? ? finish : start
    end

    protected

    def cache_extended_borders
      return if @extended_start || @extended_finish
      data.extend_borders! contig.target

      @extended_start = data.start contig.target
      @extended_finish = data.finish contig.target

      if @extended_start < 0 || @extended_finish < 0
        raise 'Wrong extension in ' + self.inspect
      end

      data.assign_attr_by_target(:start, start, contig.target)
      data.assign_attr_by_target(:finish, finish, contig.target)
    end
  end
end
