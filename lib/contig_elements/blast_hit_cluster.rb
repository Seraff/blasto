module ContigElements
  class BlastHitCluster < Basic
    FINAL_CLUSTERS_FILENAME = 'final_clusters.gff'

    def initialize(contig, start, finish, data, extra_data: {})
      super
      @start, @finish = forward? ? [start, finish].sort : [start, finish].sort.reverse
    end

    def blast_hits
      data
    end

    def blast_hits=(val)
      self.data = val
    end

    def best_blast_hit
      @best_blast_hit ||= blast_hits.first
    end

    def forward?
      [1,2,3].include? extra_data[:frame]
    end

    def right_border
      [start, finish].sort.last
    end

    def left_border
      [start, finish].sort.first
    end

    def to_gff
      notes = "ID=#{SecureRandom.hex}"
      [contig.title, :blast, :gene, start, finish, '.', '+', 0, notes].join("\t")
    end
  end
end
