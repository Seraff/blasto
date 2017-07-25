class Filterer
end

Dir["#{ROOT_PATH}/lib/helpers/filterers/*.rb"].each {|file| require file }
