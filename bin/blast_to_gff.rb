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
  o.string '-in', '(required) Input file.'
  o.string '-out', 'Output GFF file.'
  o.string '-t', '--target', '(required) Source for generating gff file (query|subject)'
  o.bool '-b', '--back_translate', 'Back translate BLAST hit coordinates'
  o.string '-m', '--mode', "Mode for input file parsing (genome|transcriptome, do nothig by default)\n #{mode_help}"
  o.bool '--extend', 'Make an extended blast hits gff file'
  o.bool '--show_extended', 'Show extended regions in .gff'
  o.bool '--merge', 'Merge close blast hits with the same query, subject and frame'
  o.integer '--max_distance', 'Max distance between close hits for merging'
  o.on '-h', '--help', 'Print options' do
    puts o
    puts "\nBy default output file is saved near the input file\n"
    puts "Examples:"
    puts "\t./blast_to_gff.rb -in hits.csv -t subject - simple conversion"
    puts "\t./blast_to_gff.rb -in hits.csv -out result.gff -t query -m transcriptome"
    puts "\t./blast_to_gff.rb -in hits.csv --target subject"
    exit
  end
end

assure_params_provided params, :in, :t

unless %w(query subject).include? params[:t]
  puts "Wrong source option, use --help for more info"
  exit 0
end

input_file_path = params[:in]
output_file_path = params[:out] || change_path(params[:in], new_ext: 'gff')

reader = BlastReader.new input_file_path

if params[:merge]
  merged_path = change_path params[:in], append: "merged_#{SecureRandom.hex}"
  merged_file = File.open merged_path, 'w'

  pb = ProgressBar.create(title: 'Merging', starting_at: 0, total: reader.hits_count)
  reader.merge_hits merged_file, target: params[:t], max_distance: params[:max_distance] || 256, progress_bar: pb
  reader.close

  reader = BlastReader.new merged_path
end

gff_file = File.open output_file_path, 'w'
gff_file.puts "##gff-version 3"

title = 'Converting'
title += ' and back-translating' if params[:b]
pb = ProgressBar.create(title: title, starting_at: 0, total: reader.hits_count)

reader.each_hit do |hit|
  hit.back_translate_coords! params[:target] if params[:b]
  hit.extend_borders! params[:target] if params[:extend]

  gff = hit.to_gff params[:target], extra_data_keys: [:evalue], show_extended: params[:show_extended]

  # additional gff processing
  gff_array = gff.split("\t")

  if params[:mode]
    case params[:mode].to_sym
    when :genome
      gff_array[0].gsub!(/_\d+\z/, '')
    when :transcriptome
      gff_array[0].gsub!(/_length_.+/, '')
    end
  end

  gff_file.puts gff_array.join("\t")

  pb.increment
end

puts "Finished"

gff_file.close
reader.file.close
File.delete(merged_file) if merged_file
