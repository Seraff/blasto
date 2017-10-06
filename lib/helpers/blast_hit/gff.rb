class BlastHit
  module Gff
    COLOR_MERGED = '#ea2502'

    def to_gff(target, extra_data_keys: [], show_extended: false, show_merged: false)
      keys = detect_keys target

      required_fields = [keys[:id], keys[:start], keys[:finish]]
      required_fields.each do |f|
        raise 'Cannot convert: not enough fields' unless data.key?(f)
      end

      result = []
      result << gff_array(target, extra_data_keys)
      result += gff_exons(target, show_extended, show_merged)

      result.map { |r| r.join("\t") }.join("\n")
    end

    def gff_exons(target, show_extended, show_merged)
      return [] if !show_extended && !show_merged

      keys = detect_keys target
      result = []
      coords = []

      exon_start, exon_finish = data[keys[:start]], data[keys[:finish]]

      if extended? && show_extended
        exon_start, exon_finish = unextended_data[keys[:start]], unextended_data[keys[:finish]]
        coords = [[exon_start, exon_finish]]
      end

      if merged? && show_merged
        coords = [[exon_start]]
        merging_gaps.each do |s, f|
          coords[-1] << s
          coords << [f]
        end
        coords[-1] << exon_finish
      end

      if extended? && show_extended
        [data[keys[:start]], data[keys[:finish]]].each do |c|
          coords << [c, c]
        end
      end

      coords.each do |s, f|
        new_id = SecureRandom.hex
        new_note = "ID=#{new_id};Parent=#{gff_id(target)}"
        s, f = [s, f].sort
        result << [data[keys[:id]], 'blast', 'exon', s, f, '.', gff_strand(target), gff_frame(target), new_note]
      end

      result
    end

    def gff_array(target, extra_data_keys = [])
      keys = detect_keys target
      contig_name = data[keys[:id]]
      s, f = [data[keys[:start]], data[keys[:finish]]].sort
      s = 1 if s.zero?
      [contig_name, 'blast', 'gene', s, f, '.', gff_strand(target), gff_frame(target), gff_notes(target, extra_data_keys)]
    end

    def gff_notes(target, extra_data_keys = [])
      keys = detect_keys target
      notes = "ID=#{gff_id(target)}"

      if extra_data_keys.any?
        extra_data = extra_data_keys.select { |k| data[k] }.map { |k| "#{k}=#{data[k]}" }.join(', ')
        notes += ";Note=#{extra_data}"
      end

      notes += ";Evalue=#{data[:evalue]}"
      notes += ";Color=#{color}" if color
      notes
    end

    def gff_frame(target)
      detect_frame target
    end

    def gff_strand(target)
      keys = detect_keys target
      data[keys[:start]] < data[keys[:finish]] ? '+' : '-'
    end

    def gff_id(target)
      keys = detect_keys target

      @gff_id ||= {}
      @gff_id[data[keys[:opposite_id]]] ||= begin
        [data[keys[:opposite_id]], (1..16).to_a.map { (0..9).to_a.sample }.join, gff_frame(target)].join('_')
      end
    end

    def color
      merged? ? COLOR_MERGED : nil
    end
  end
end
