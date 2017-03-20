require 'zenity'

module Helpers
  CODON_TABLE = {"TTT"=>"F", "TTC"=>"F", "TTA"=>"L",
                 "TTG"=>"L", "TCT"=>"S", "TCC"=>"S",
                 "TCA"=>"S", "TCG"=>"S", "TAT"=>"Y",
                 "TAC"=>"Y", "TAA"=>"Stop", "TAG"=>"Stop",
                 "TGT"=>"C", "TGC"=>"C", "TGA"=>"Stop",
                 "TGG"=>"W", "CTT"=>"L", "CTC"=>"L",
                 "CTA"=>"L", "CTG"=>"L", "CCT"=>"P",
                 "CCC"=>"P", "CCA"=>"P", "CCG"=>"P",
                 "CAT"=>"H", "CAC"=>"H", "CAA"=>"Q",
                 "CAG"=>"Q", "CGT"=>"R", "CGC"=>"R",
                 "CGA"=>"R", "CGG"=>"R", "ATT"=>"I",
                 "ATC"=>"I", "ATA"=>"I", "ATG"=>"M",
                 "ACT"=>"T", "ACC"=>"T", "ACA"=>"T",
                 "ACG"=>"T", "AAT"=>"N", "AAC"=>"N",
                 "AAA"=>"K", "AAG"=>"K", "AGT"=>"S",
                 "AGC"=>"S", "AGA"=>"R", "AGG"=>"R",
                 "GTT"=>"V", "GTC"=>"V", "GTA"=>"V",
                 "GTG"=>"V", "GCT"=>"A", "GCC"=>"A",
                 "GCA"=>"A", "GCG"=>"A", "GAT"=>"D",
                 "GAC"=>"D", "GAA"=>"E", "GAG"=>"E",
                 "GGT"=>"G", "GGC"=>"G", "GGA"=>"G",
                 "GGG"=>"G"}

  def select_file(subject: nil, title: nil)
    unless title
      subject = subject ? ": #{subject}" : subject
      title = "Choose a file#{subject}"
    end
    a = Zenity::send(:'file-selection', title: title)
    a.strip
  end

  def open_file(subject: nil, title: nil)
    filename = select_file(subject: subject, title: title)
    raise 'You should specify a file!' unless File.exists?(filename)
    File.open filename, 'r'
  end

  def append_to_filename(filename, what)
    dir, file = Pathname.new(filename).split
    file = file.to_s.split('.')
    file[0] += "_#{what}"
    file = file.join('.')
    (dir + Pathname.new(file)).to_s
  end

  def change_path(path, new_dir:, append:, new_ext: nil)
    base_name = Pathname.new(path).basename.to_s
    base_name = append_to_filename(base_name, append)
    base_name = base_name.split('.')[0] + ".#{new_ext}" if new_ext
    (Pathname.new(new_dir) + Pathname.new(base_name)).to_s
  end

  def assure_params_provided(params, *options)
    options.each do |o|
      if !params[o]
        puts "Not all required options provided, use --help for more info"
        exit
      end
    end
  end
end

include Helpers

class Bio::GFF::GFF3::Record
  def reverse_strand?
    raise "Frame is not provided: node #{seqname}" unless frame
    [4, 5, 6].include? frame
  end
end
