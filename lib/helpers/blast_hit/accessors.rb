class BlastHit
  module Accessors
    ACCESSOR_METHODS = [:id, :start, :finish, :frame, :len]

    ACCESSOR_METHODS.each do |attr_name|
      define_method attr_name do |target = nil|
        target = target || @target_context
        raise 'Target is not specified' unless target

        data[TARGET_KEYS[target.to_sym][attr_name.to_sym]]
      end

      define_method :"#{attr_name}=" do |val|
        target = @target_context
        raise 'Target is not specified' unless target

        @data[TARGET_KEYS[target.to_sym][attr_name.to_sym]] = val
      end
    end

    def with_target_context(target)
      @target_context = target
      yield self
      @target_context = nil
    end
  end
end
