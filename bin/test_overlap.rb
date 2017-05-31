#!/usr/bin/env ruby
require_relative 'lib.rb'

class Element
  attr_accessor :start, :finish

  def initialize(start, finish)
    @start = start
    @finish = finish
  end

  def intersects?(other)
    IntervalsHelper.intersects? start..finish, other.start..other.finish
  end
end

def hit_clusters
  [Element.new(50, 100), Element.new(10, 30), Element.new(25, 55), Element.new(120, 150)]
end

def group_blasters_by_intersections
  groups = []
  sorted = hit_clusters.dup.sort_by { |h| h.start }

  current_group = []
  sorted.each do |e|
    found = false

    if current_group.empty?
      current_group << e
      next
    end

    current_group.each do |ge|
      if e.intersects?(ge)
        current_group << e
        found = true
        break
      end
    end

    unless found
      groups << current_group
      current_group = [e]
    end
  end

  groups << current_group

  groups
end

puts group_blasters_by_intersections.inspect
