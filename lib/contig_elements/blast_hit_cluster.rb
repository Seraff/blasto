module ContigElements
  class BlastHitCluster < Basic
    FINAL_CLUSTERS_FILENAME = 'final_clusters.gff'

    def initialize(contig, start, finish, data, extra_data: {})
      super
      if [1,2,3].include? extra_data[:frame]
        @start, @finish = [start, finish].sort
      else
        @start, @finish = [start, finish].sort.reverse
      end
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

    def to_gff
      notes = "ID=#{SecureRandom.hex}"
      [contig.title, :blast, :gene, start, finish, '.', '+', 0, notes].join("\t")
    end
  end
end
