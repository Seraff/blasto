#!/usr/bin/env ruby
require_relative '../lib.rb'

fasta_string = <<END_OF_STRING
>gi|398365175|ref|NP_009718.3| Cdc28p [Saccharomyces cerevisiae S288c]
MSGELANYKRLEKVGEGTYGVVYKALDLRPGQGQRVVALKKIRLESEDEGVPSTAIREISLLKELKDDNI
VRLYDIVHSDAHKLYLVFEFLDLDLKRYMEGIPKDQPLGADIVKKFMMQLCKGIAYCHSHRILHRDLKPQ
NLLINKDGNLKLGDFGLARAFGVPLRAYTHEIVTLWYRAPEVLLGGKQYSTGVDTWSIGCIFAEMCNRKP
IFSGDSEIDQIFKIFRVLGTPNEAIWPDIVYLPDFKPSFPQWRRKDLSQVVPSLDPRGIDLLDKLLAYDP
INRISARRAAIHPYFQES
END_OF_STRING

seq = Bio::FastaFormat.new fasta_string
contig = Contig.new(seq)

## SL
contig.sl_mappings = ContigElementCollection.new

contig.sl_mappings << ContigElements::Basic.new(contig, 60, 61, nil)
contig.sl_mappings << ContigElements::Basic.new(contig, 135, 136, nil)

coll = ContigElementCollections::Zoi.new
zoi = ContigElements::Zoi.new(contig, 50, 150, '')
coll << zoi

puts "#{zoi.start} - #{zoi.finish}"

coll.shorten_to_sl

puts "#{zoi.start} - #{zoi.finish}"

