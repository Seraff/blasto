require 'bio'
require 'helpers'

class Annotator
  attr_accessor :genome

  @@hits_prepared = false

  def self.prepare_blast_hits
    return true if @@hits_prepared

    ## take blast_hit file
    ## split by contigs
    ## save each contig to it's own file

    ##
    @@hits_prepared = true
  end

  def initialize
    file_name = select_file('genome contigs')
    @genome = Bio::FlatFile.open(Bio::FastaFormat, file_name)
  end

  def each_contig
    @genome.each do |fasta_format|
      yield Contig.new(fasta_format)
    end
  end

  def annotate
    each_contig { |c| c.annotate }
  end
end
