class DummyContig
  attr_accessor :sl_mappings, :zoi, :blast_hit_clusters, :seq

  def initialize(sl_data, zoi_data, hit_clusters_data, seq:)
    self.seq = seq || _generate_seq

    sls = sl_data.map { |e| ContigElements::Sl.new self, e[0], e[1], { coverage: e[2] } }
    self.sl_mappings = ContigElementCollection.new(sls)

    elements = zoi_data.map do |z|
      ContigElements::Zoi.new(self, z[0], z[1], '')
    end

    self.zoi = ContigElementCollections::Zoi.new(elements)
    self.zoi.contig = self if self.zoi.any?

    clusters = hit_clusters_data.map do |e|
      extra_data = { frame: e[2], forward: [1, 2, 3].include?(e[2]) }
      ContigElements::BlastHitCluster.new(self, e[0], e[1], [], extra_data: extra_data)
    end

    self.blast_hit_clusters = ContigElementCollection.new(clusters)
  end

  def title
    'test_contig'
  end

  def length
    seq.length
  end

  protected

  def _generate_seq
    (0..1000).map { |e| %w(G, T, A, C).sample }.join
  end
end

class TestDataset
  attr_accessor :built, :contig, :zoi_coords,
    :sl_coords, :hit_cluster_coords, :seq

  def initialize(contig_seq = nil)
    self.seq = contig_seq
    self.contig = nil
    self.zoi_coords = []
    self.sl_coords = []
    self.hit_cluster_coords = []
    self.built = false
  end

  def add_zoi(*coords)
    self.zoi_coords += coords
    self.built = false
  end

  def add_sl(*coords)
    self.sl_coords += coords
    self.built = false
  end

  def add_hit_cluster(*coords)
    self.hit_cluster_coords += coords
    self.built = false
  end

  def contig
    unless built
      rebuild
      built = true
    end

    @contig
  end

  protected

  def rebuild()
    self.contig = DummyContig.new(self.sl_coords,
                                  self.zoi_coords,
                                  self.hit_cluster_coords,
                                  seq: seq)
  end
end
