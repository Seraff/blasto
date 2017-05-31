module IntervalsHelper
	class << self
		def intersects?(a, b)
	    ![a, b].intersection.nil?
		end

		def covers?(a, b)
	    [a, b].intersection == (b.first..b.last)
		end
	end
end
