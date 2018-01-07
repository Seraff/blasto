module ContigElements
  class Sl < Basic
    def coverage
      @coverage ||= data[:coverage].to_f
    end

    def center_coord_for_frame(frame)
      return (start+finish)/2 if na_len > 2

      if [1, 2, 3].include?(frame)
        finish
      elsif [4, 5, 6].include?(frame)
        start
      else
        raise "Wrong frame for ContigSubsequence#forward_frame?(): #{frame}"
      end
    end
  end
end
