#!/usr/bin/env ruby
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
  o.bool '--merge', 'Merge close blast hits with the same query, subject and frame'
  o.bool '--extend', 'Make an extended blast hits gff file'
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

if params[:merge]
  merged_path = change_path(params[:in], append: "merged_#{SecureRandom.hex}")
  merged_file = File.open(merged_path, 'w')

  pb = ProgressBar.create(title: 'Merging', starting_at: 0, total: reader.hits_count)
  reader.merge_hits merged_file, target: params[:t], max_distance: params[:max_distance] || 256, progress_bar: pb

  reader.close

  reader = BlastReader.new merged_path
end

outfile_path = params[:out]
outfile_path ||= append_to_filename(params[:in], 'non_translated').gsub(/\.\w+\z/, '.gff')

gff_file = File.open outfile_path, 'w'
gff_file.puts "##gff-version 3"

pb = ProgressBar.create(title: 'Translating hits', starting_at: 0, total: reader.hits_count)

reader.back_translate_to_gff gff_file, target: params[:target],
                                       mode: params[:mode],
                                       progress_bar: pb,
                                       extend_borders: params[:extend]

puts "Finished"

gff_file.close
reader.file.close
File.delete(merged_file) if merged_file
