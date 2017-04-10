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

  def change_path(path, new_dir: nil, append: nil, new_ext: nil)
    pathname = Pathname.new(path)
    base_name = pathname.basename.to_s
    new_dir ||= pathname.dirname.to_s

    base_name = append_to_filename(base_name, append) if append
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

  def assure_file_param_has_extension(params, key, extension)
    required_ext = extension.to_s.gsub('.', '')
    actual_ext = Pathname.new(params[key.to_sym]).extname.gsub('.', '')

    return true if required_ext == actual_ext

    puts "File from option #{key} has incorrect extension"
    exit
  end
end

include Helpers

class Bio::GFF::GFF3::Record
  def reverse_strand?
    raise "Frame is not provided: node #{seqname}" unless frame
    [4, 5, 6].include? frame
  end
end

class Bio::Sequence::NA
  def aa_codon_usage
    usage = {}

    codon_usage.each do |codon, count|
      aa = Bio::Sequence::NA.new(codon).translate
      usage[aa] ||= {}
      usage[aa][codon.upcase] ||= 0
      usage[aa][codon.upcase] += count
    end

    usage
  end

  def aa_codon_usage_statistics
    usage = aa_codon_usage
    new_usage = {}

    aa_count = usage.values.map(&:values).flatten.inject(0, :+)

    usage.each do |aa_name, codons|
      codons_count = codons.values.flatten.inject(0, :+)

      new_usage[aa_name] ||= {}
      new_usage[aa_name] = codons.map do |codon, count|
        stats = { count: count,
                  total_percentage: count.to_f/aa_count,
                  aa_percentage: count.to_f/codons_count }
        [codon, stats]
      end.to_h
    end

    new_usage
  end
end
