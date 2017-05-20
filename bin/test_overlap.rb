#!/usr/bin/env ruby
require_relative 'lib.rb'

def intersects?(a, b)
  ![a, b].intersection.nil?
end

def covers?(a, b)
  [a, b].intersection == (b.first..b.last)
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

def prepare_intervals(data)
  results = { good: [], bad: [] }

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

# a = [
#   221_542..226_510,
#   218_232..221_627,
#   226_728..227_342,
#   226_476..226_805,
#   227_955..229_299
# ]

a = [
  500..1000,
  100..600,
  1500..2000,
  800..1800,
  2500..3000,
  3500..4000,
  3800..4500,
  2600..2900
]

a = [221_542..226_510, 218_232..221_627, 226_728..227_342, 226_476..226_805, 227_955..229_299]
filter_covered a
r = prepare_intervals a
puts r.inspect
