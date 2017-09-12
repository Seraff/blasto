module Filterers
  class BestNonoverlappedBlastHits < Filterer
    def initialize(input_path, output_path: nil, target: nil)
      @reader = BlastReader.new input_path
      @output_path = output_path || input_path
      @target = target
    end

    def filtered?
      !@output_hits.nil?
    end

    def perform
      filter_dict
      save_result
    end

    def filter_dict
      heap = @reader.hits.dup
      bests = []

      while heap.any? do
        best = get_best(heap)
        intersected = get_intersected(heap, best)

        heap -= (intersected + [best])
        bests << best
      end

      @output_hits = bests
    end

    def save_result
      @reader.modify_hits(@output_path) do |h|
        @output_hits.include? h
      end

      @reader.close
    end

    protected

    def intersects_with?(first, second)
      begin
        first = [first.attr_by_target(:start, @target), first.attr_by_target(:finish, @target)].sort
      rescue
        puts [first.attr_by_target(:start, @target), first.attr_by_target(:finish, @target)].inspect
        exit
      end
      second = [second.attr_by_target(:start, @target), second.attr_by_target(:finish, @target)].sort
      IntervalsHelper.intersects? first[0]..first[1], second[0]..second[1]
    end

    def get_intersected(els, element)
      els.select { |e| e != element && intersects_with?(e, element) }
    end

    def get_best(els)
      els.max { |a, b| [b.data[:evalue].to_f, a.data[:bitscore].to_f] <=> [a.data[:evalue].to_f, b.data[:bitscore].to_f] }
    end
  end
end
