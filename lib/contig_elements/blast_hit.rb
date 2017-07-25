module ContigElements
  class BlastHit < Basic
    def frame
      @frame ||= data.frame contig.target
    end

    def organism
      @organism ||= data.data[:qseqid].split(/_|\|/).first
    end

    def gene
      @gene ||= data.data[:qseqid].split(/_|\|/)[1..-1].join
    end

    def extended_start
      cache_extended_borders
      @extended_start
    end

    def extended_finish
      cache_extended_borders
      @extended_finish
    end

    protected

    def cache_extended_borders
      return if @extended_start && @extended_finish
      data.extend_borders! contig.target

      @extended_start = data.start contig.target
      @extended_finish = data.finish contig.target

      data.assign_attr_by_target(:start, start, contig.target)
      data.assign_attr_by_target(:finish, finish, contig.target)
    end
  end
end
