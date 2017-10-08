## Hardcoded system paths

class Preparer
  module Paths
    SPLITTED_DATA_FOLDER = 'contigs'

    PATHS = {
      hits_csv: 'hits.csv',
      hits_gff: 'hits.gff',
      transcripts_csv: 'transcripts.csv',
      transcripts_gff: 'transcripts.gff',
      clusters_csv: 'hit_clusters.csv',
      clusters_gff: 'hit_clusters.gff',
      clusters_extended_gff: 'hit_clusters_extended.gff',
      transcripts_bin: 'transcripts_bin.csv',
      sl_mapping: 'sl_mapping.bed',
      back_translated_hits: 'hits_back_translated.csv',
      after_merging_hits: 'hits_after_merging.csv',
      merged_hits_gff: 'hits_merged.gff',
      merged_nonoverlapped_hits_gff: 'hits_merged_nonoverlapped.gff',
      extended_hits: 'hits_extended.csv'
    }

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

      PATHS.each do |name, filename|
        define_method("#{name}_path") do |folder|
          contig_folder_path folder, filename: filename
        end
      end

      def sorted_by_target_path(what)
        abs_path_for "#{what}_sorted_by_target.csv"
      end
    end
  end
end
