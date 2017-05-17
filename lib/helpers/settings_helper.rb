class SettingsHelper
  @@instance = SettingsHelper.new

  def self.instance
    @@instance
  end

  def tmp_abs_pathname
    abs_path_for_setting :tmp_dir
  end

  def abs_path_for_setting(name)
    make_abs_ath Pathname.new(Settings.annotator.send(name))
  end

  protected

  def make_abs_ath(rel_path)
    rel_path.absolute? ? rel_path : Pathname.new(ROOT_PATH).join(rel_path)
  end
end
