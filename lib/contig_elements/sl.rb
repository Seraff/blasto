module ContigElements
  class Sl < Basic
    def coverage
      @coverage ||= data[:coverage].to_f
    end
  end
end
