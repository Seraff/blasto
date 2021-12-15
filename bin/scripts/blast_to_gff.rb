#!/usr/bin/env ruby
require_relative 'lib.rb'

mode_help = %q{
        result seqid names are formed by the rules:
        spades mode: SEQID_length_<nucleotides count>[_optional_string]
        short mode:  SEQID

        Example:
        seqid:                NODE_1096_length_1211_cov_913.715_1
        cound of nucleotides: 1211
        frame:                1

        gff seqid (spades mode): NODE_1096_length_1211_cov_913.715
        gff seqid (short mode):  NODE_1096
}

banner = %q{
  Convert blast hits file with amino acid coordinates to gff annotation with nucleotide indicies

  seqid must contain nt length and frame: SEQID_length_<nucleotides count>[_optional_string]_<frame>
}


params = Slop.parse do |o|
  o.banner = ""
  o.string '-in', '(required) Input file'
  o.string '-out', 'Output GFF file (same name and folder as input by default)'
  o.string '-t', '--target', '(required) Source for generating gff file (query|subject)'
  o.bool '-b', '--back_translate', 'Back translate BLAST hit coordinates'
  o.string '-m', '--mode', "Strategy for ouput seqid formatting (spades|short, do not change seqid by default)\n #{mode_help}"
  o.bool '--extend', 'Make an extended blast hits gff file'
  o.bool '--show_extended', 'Show extended regions in .gff (using intron/exon notation). Ignored if --extend option not provided.'
  o.bool '--merge', 'Merge close blast hits with the same query, subject and frame'
  o.bool '--show_merged', 'Show merged hits in .gff (using intron/exon notation). Ignored if --extend option not provided.'
  o.integer '--max_distance', 'Max distance between close hits for merging (in na or aa, depending on the source file). Ignored if --merge option not provided. '
  o.on '-h', '--help', 'Print options' do
    puts o
    puts "\nBy default output file is saved near the input file\n"
    puts "Examples:"
    puts "\t./blast_to_gff.rb -in hits.csv -t subject - simple conversion"
    puts "\t./blast_to_gff.rb -in hits.csv -out result.gff -t query -m short"
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
reader.cache_hits

if params[:merge]
  pb = ProgressBar.create(title: 'Merging', starting_at: 0, total: reader.hits_count)
  default_max_distance = params[:back_translate] ? 85 : 256
  reader.merge_hits! target: params[:t], max_distance: params[:max_distance] || default_max_distance, progress_bar: pb
end

if params[:extend]
  pb = ProgressBar.create(title: 'Extending', starting_at: 0, total: reader.hits_count)
  reader.extend_borders! target: params[:target], progress_bar: pb
end

if params[:back_translate]
  pb = ProgressBar.create(title: 'Back translating', starting_at: 0, total: reader.hits_count)
  reader.each_hit do |hit|
    hit.back_translate_coords! params[:target]
    pb.increment
  end
end

File.open(output_file_path, 'w') do |f|
  pb = ProgressBar.create(title: 'Writing to file', starting_at: 0, total: reader.hits_count)
  f.puts "##gff-version 3"

  reader.each_hit do |hit|
    show_extended = params[:extend] ? params[:show_extended] : false
    show_merged = params[:merge] ? params[:show_merged] : false

    gff = hit.to_gff params[:target], extra_data_keys: [:evalue], show_extended: show_extended, show_merged: show_merged
    gff_array = gff.split("\t")

    if params[:mode]
      case params[:mode].to_sym
      when :spades
        gff_array[0].gsub!(/_\d+\z/, '')
      when :short
        gff_array[0].gsub!(/_length_.+/, '')
      end
    end

    f.puts gff_array.join("\t")
    pb.increment
  end
end

reader.close
