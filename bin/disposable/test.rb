#!/usr/bin/env ruby
require_relative '../lib.rb'
require 'matplotlib/pyplot'

plt = Matplotlib::Pyplot

xs = [*1..100].map {|x| (x - 50) * Math::PI / 100.0 }
ys = xs.map {|x| Math.sin(x) }

xs = [*-10..10]
ys = xs.map { |x| x*2 }

plt.plot(xs, ys)
plt.show()
