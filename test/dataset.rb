class DummyContig
  attr_accessor :sl_mappings, :zoi, :blast_hits, :seq

  def initialize(sl_data, zoi_data, blast_hits_data, seq:)
    self.seq = seq || _generate_seq

    sls = sl_data.map { |e| ContigElements::Sl.new self, e[0], e[1], { coverage: e[2] } }
    self.sl_mappings = ContigElementCollection.new(sls)

    elements = zoi_data.map do |z|
      ContigElements::Zoi.new(self, z[0], z[1], '')
    end

    self.zoi = ContigElementCollections::Zoi.new(elements)
    self.zoi.contig = self if self.zoi.any?

    hits = blast_hits_data.map do |e|
      data = { qseqid: 'guy|12345', sframe: e[2] }
      hit = BlastHit.new(data.keys, data)
      ContigElements::BlastHit.new(self, e[0], e[1], hit)
    end

    self.blast_hits = ContigElementCollection.new(hits)
  end

  def title
    'test_contig'
  end

  def length
    seq.length
  end

  def target
    :subject
  end

  def subsequence(left, right)
    ContigSubsequence.new self, left, right
  end

  protected

  def _generate_seq
    (0..1000).map { |e| %w(G, T, A, C).sample }.join
  end
end

class TestDataset
  attr_accessor :built, :contig, :zoi_coords,
    :sl_coords, :blast_hit_coords, :seq

  def initialize(contig_seq = nil)
    self.seq = contig_seq
    self.contig = nil
    self.zoi_coords = []
    self.sl_coords = []
    self.blast_hit_coords = []
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

  def add_hit(*coords)
    self.blast_hit_coords += coords
    self.built = false
  end

  def contig
    unless built
      rebuild
      self.built = true
    end

    @contig
  end

  protected

  def rebuild()
    self.contig = DummyContig.new(self.sl_coords,
                                  self.zoi_coords,
                                  self.blast_hit_coords,
                                  seq: seq)
  end
end
