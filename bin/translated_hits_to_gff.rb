#!/usr/bin/ruby
require_relative 'lib.rb'
mode_help = %q{
        qseqid/sseqid must be in format: SEQID_length_<nucleotides count>[_optional_string]_<frame>
        result seqid names are formed by the rules
        genome mode (spades ids): SEQID_length_<nucleotides count>[_optional_string]
        transcriptome mode (other ids) : SEQID

        Example:
        seqid: NODE_1096_length_1211_cov_913.715_1
        cound of nucleotides: 1211
        frame: 1
        gff seqid (genome mode): NODE_1096_length_1211_cov_913.715
        gff seqid (transcriptome mode): NODE_1096
}

params = Slop.parse do |o|
  o.banner = "Convert blast hits file with amino acid coordinates to gff annotation with nucleotide indicies\n"
  o.string '-in', '--input_file', 'Input file. If not provided, opens the file selection window.'
  o.string '-out', '--output_file', 'Output file. If not provided, puts the output gff file to the input file folder.'
  o.string '-t', '--target', '(required) Source for generating gff file (query|subject)'
  o.string '-m', '--mode', "(required) Mode for input file parsing (genome|transcriptome)\n #{mode_help}"
  o.on '-h', '--help', 'Print options' do
    puts o
    puts "\nBy default output file is saved near the input file\n"
    puts "Examples:"
    puts "\t./translated_hits_to_gff.rb -in hits.csv -out result.gff -t query -m genome"
    puts "\t./translated_hits_to_gff.rb -in hits.csv --target subject --mode transcriptome"
    puts "\t./translated_hits_to_gff.rb -t query -m genome"
    exit
  end
end

[:t, :m].each do |opt|
  unless params[opt]
    puts "-#{opt} option is required, use --help for more info"
    exit
  end
end

input_file_path = params[:in] || select_file(subject: "translated blast hits")
reader = BlastReader.new input_file_path

outfile_path = params[:out]
outfile_path ||= append_to_filename(reader.file.path, 'non_tranlated').gsub(/\.\w+\z/, '.gff')

gff_file = File.open outfile_path, 'w'
gff_file.puts "##gff-version 3"

pb = ProgressBar.create(title: 'Translating hits', starting_at: 0, total: reader.hits_count)

reader.each_hit do |hit|
  hit.back_translate_coords! params[:target]
  gff = hit.to_gff params[:target]

  # additional gff processing
  gff_array = gff.split("\t")

  case params[:mode].to_sym
  when :genome
    gff_array[0].gsub!(/_\d+\z/, '')
  when :transcriptome
    gff_array[0].gsub!(/_length_.+/, '')
  end

  gff_file.puts gff_array.join("\t")

  pb.increment
end

puts "Finished"

gff_file.close
reader.file.close
