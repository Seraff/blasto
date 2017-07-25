module ContigElements
  class BlastHitCluster < Basic
    FINAL_CLUSTERS_FILENAME = 'final_clusters.gff'

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
