## IMPORTANT: for test purposes only

#!/usr/bin/ruby
require_relative 'lib.rb'

# inputs
genome_path = 'data/test/p57_DNA_nucleotides.fa'
hits_path = 'data/Trinity-GG_p57_6_frames_translat_bl_report_best.csv'
reads_path = 'data/datasets/p57_RNA_to_DNA/p57_bw2_sorted.bam'

annotator = Annotator.new genome_path: genome_path,
                          hits_path: hits_path,
                          reads_path: reads_path
