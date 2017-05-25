require_relative 'blast_reader'

class TranscriptomeCleaner
  attr_accessor :data_hash, :target, :reader, :folder
  BIN_PATH = 'tmp/annotator/transcriptome_bin.csv'

  def initialize(folder, target:)
    @folder = folder
    @target = target
    @reader = BlastReader.new Preparer.transcripts_csv_path(folder)
    @reader.cache_hits
  end

  def clean
    normalize

    filter_by_size
    filter_by_intersections
  end

  protected

  def normalize
    # normalizing hits: 2..1 => 1..2
    @reader.each_hit do |hit|
      hit.with_target_context(target) do |h|
        if h.start > h.finish
          old_start = h.start
          h.start = h.finish
          h.finish = old_start
        end
      end
    end
  end

  def filter_by_size
    short_hits = []

    @reader.hits.keep_if do |hit|
      result = (hit.finish(target) - hit.start(target) + 1) >= Settings.annotator.transcriptome_min_size
      short_hits << hit unless result
      result
    end

    BadTranscriptsLogger.add_hit_collection_to_bin(folder, hits: short_hits, reason: :short)
  end

  def filter_by_intersections
    data_hash = @reader.hits.map do |h|
      [(h.start(target))..(h.finish(target)), h]
    end
    data_hash = data_hash.to_h

    result = classify_intervals data_hash.keys
    result = result.map { |k, v| [k, v.map { |r| data_hash[r] }] }.to_h

    BadTranscriptsLogger.add_hit_collection_to_bin(folder, hits: result[:bad], reason: :intersected)

    result
  end

  def classify_intervals(data)
    data = data.uniq
    results = { good: [], bad: [] }

    # clearing covered intervals
    data = filter_covered data

    # detecting absolutely good groups
    data.each do |me|
      next if results[:bad].include? me

      without_me = data.dup
      without_me.delete(me)
      intersection_found = false

      without_me.each do |other|
        if intersects? other, me
          intersection_found = true
          results[:bad] += [me, other]
          break
        end
      end

      unless intersection_found
        results[:good] << me
      end
    end

    results
  end

  def filter_covered(data)
    covered = []

    data.each do |i|
      add_to_covered = false

      other = data.dup
      other.delete(i)

      other.each do |j|
        if covers? j, i
          add_to_covered = true
          break
        end
      end

      covered << i if add_to_covered
    end

    data - covered
  end

  def intersects?(a, b)
    ![a, b].intersection.nil?
  end

  def covers?(a, b)
    [a, b].intersection == (b.first..b.last)
  end
end
