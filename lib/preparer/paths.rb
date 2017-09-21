## Hardcoded system paths

class Preparer
  module Paths
    SPLITTED_DATA_FOLDER = 'contigs'
    # per-folder files
    CSV_HITS_FILENAME = 'hits.csv'
    GFF_HITS_FILENAME = 'hits.gff'
    CSV_TRANSCRIPTS_FILENAME = 'transcripts.csv'
    GFF_TRANSCRIPTS_FILENAME = 'transcripts.gff'
    CSV_CLUSTERS_FILENAME = 'hit_clusters.csv'
    GFF_CLUSTERS_FILENAME = 'hit_clusters.gff'
    GFF_CLUSTERS_EXTENDED_FILENAME = 'hit_clusters_extended.gff'
    TRANSCRIPTS_BIN_FILENAME = 'transcripts_bin.csv'
    SL_MAPPING_FILENAME = 'sl_mapping.bed'

    # global files
    BACK_TRANSLATED_HITS_FILENAME = 'hits_back_translated.csv'
    MERGED_HITS_FILENAME = 'hits_merged.csv'
    EXTENDED_HITS_FILENAME = 'hits_extended.csv'

    def self.included base
      base.extend ClassMethods
    end

    module ClassMethods
      def tmp_abs_pathname
        SettingsHelper.instance.tmp_abs_pathname
      end

      def abs_path_for(name)
        tmp_abs_pathname + Pathname.new(name)
      end

      def contigs_folder_path
        abs_path_for SPLITTED_DATA_FOLDER
      end

      def contig_folder_path(folder, filename: nil)
        path = contigs_folder_path + Pathname.new(folder)
        path += Pathname.new(filename) if filename
        path
      end

      def sorted_by_target_path(what)
        abs_path_for "#{what}_sorted_by_target.csv"
      end

      def back_translated_hits_path(folder)
        contig_folder_path folder, filename: BACK_TRANSLATED_HITS_FILENAME
      end

      def merged_hits_path(folder)
        contig_folder_path folder, filename: MERGED_HITS_FILENAME
      end

      def extended_hits_path(folder)
        contig_folder_path folder, filename: EXTENDED_HITS_FILENAME
      end

      def hits_csv_path(folder)
        contig_folder_path folder, filename: CSV_HITS_FILENAME
      end

      def hits_gff_path(folder)
        contig_folder_path folder, filename: GFF_HITS_FILENAME
      end

      def hit_clusters_csv_path(folder)
        contig_folder_path folder, filename: CSV_CLUSTERS_FILENAME
      end

      def hit_clusters_gff_path(folder)
        contig_folder_path folder, filename: GFF_CLUSTERS_FILENAME
      end

      def hit_clusters_extended_gff_path(folder)
        contig_folder_path folder, filename: GFF_CLUSTERS_EXTENDED_FILENAME
      end

      def transcripts_csv_path(folder)
        contig_folder_path folder, filename: CSV_TRANSCRIPTS_FILENAME
      end

      def transcripts_gff_path(folder)
        contig_folder_path folder, filename: GFF_TRANSCRIPTS_FILENAME
      end

      def transcripts_bin_path(folder)
        contig_folder_path folder, filename: TRANSCRIPTS_BIN_FILENAME
      end

      def sl_mapping_path(folder)
        contig_folder_path folder, filename: SL_MAPPING_FILENAME
      end
    end
  end
end
