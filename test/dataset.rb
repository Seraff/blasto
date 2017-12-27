class DummyContig
  attr_accessor :sl_mappings, :zoi

  def initialize(sl_data, zoi_data)
    sls = sl_data.map { |e| ContigElements::Sl.new self, e[0], e[1], { coverage: e[2] } }
    self.sl_mappings = ContigElementCollection.new(sls)

    elements = zoi_data.map do |z|
      ContigElements::Zoi.new(self, z[0], z[1], '')
    end

    self.zoi = ContigElementCollections::Zoi.new(elements)
    self.zoi.contig = self if self.zoi.any?
  end

  def blast_hit_clusters
    ContigElementCollection.new []
  end

  def title
    'test_contig'
  end

  def seq
    (0..1000).map { |e| %w(G, T, A, C).sample }.join
  end
end

class TestDataset
  attr_accessor :contig, :zoi_coords, :sl_coords

  def initialize()
    self.contig = nil
    self.zoi_coords = []
    self.sl_coords = []
  end

  def add_zoi(*coords)
    self.zoi_coords += coords
    rebuild
  end

  def add_sl(*coords)
    self.sl_coords += coords
    rebuild
  end

  protected

  def rebuild()
    self.contig = DummyContig.new(self.sl_coords, self.zoi_coords)
  end
end
