require 'hashie'

class BlastHitData < Hash
  include Hashie::Extensions::MergeInitializer
  include Hashie::Extensions::IndifferentAccess

  def initialize(hsh = {})
    super hsh.map { |k, v| [k, format_value(v)] }.to_h
  end

  protected

  def format_value(val)
    return val unless val.is_a?(String)

    if val.match(/\A\d+\z/)
      val.to_i
    elsif val.match(/\A\d+\.\d+\z/)
      val.to_f
    else
      val
    end
  end
end
