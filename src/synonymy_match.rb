#!/usr/bin/ruby
require 'optparse'
require 'fuzzy_match'
require 'amatch'
FuzzyMatch.engine = :amatch
require 'csv'
require 'set'
require 'parallel'
require 'ruby-progressbar'

options = {}
OptionParser.new do |opts|
  #Defaults
  options[:tpl] = nil
  options[:unique] = nil
  options[:gbif] = nil
  options[:gbif] = nil
  options[:threads] = 20
  options[:tpl_cache] = nil
  
  opts.banner = "tpl-grep: Match a set of input names to a TPL dump\nUsage: tpl-dwn.rb [options]"
  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

  opts.on("-t FILE", "--tpl FILE", "File with TPL synonymies (from Dylan's code)") {|x| options[:tpl] = x.to_s}
  opts.on("-u FILE", "--unique FILE", "File with unique names to match") {|x| options[:unique] = x.to_s}
  opts.on("-g FILE", "--gbif FILE", "GBIF file to be replaced") {|x| options[:gbif] = x.to_s}
  opts.on("-o FILE", "--output FILE", "Name of file to be written") {|x| options[:output] = x.to_s}
  opts.on("-n NUMBER", "--nthreads NUMBER", "Number of threads for fuzzy matching (default: 20)") {|x| options[:threads] = x.to_i}
  opts.on("-c FILE", "--cache FILE", "Where to save fuzzy matching lookup") {|x| options[:tpl_cache] = x.to_s}
  opts.on("-")
end.parse!

# Argument handling
if options[:tpl].nil? or not File.exist? options[:tpl]
  puts "TPL file missing or not specified; exiting"
  exit(false)
end
if options[:unique].nil? or not File.exist? options[:unique]
  puts "TPL file missing or not specified; exiting"
  exit(false)
end
if options[:gbif].nil? or not File.exist? options[:gbif]
  puts "TPL file missing or not specified; exiting"
  exit(false)
end
if options[:output].nil? or not File.exist? options[:output]
  puts "TPL file missing or not specified; exiting"
  exit(false)
end
if options[:tpl_cache].nil? or not File.exist? options[:tpl_cache]
  puts "TPL cache file missing or not specified; exiting"
  exit(false)
end

# Load synonyms and build (fuzzy) match table
tpl_lookup = Hash.new
CSV.foreach(options[:tpl]) do |line|
  line.each {|sp| tpl_lookup[sp.downcase.sub("_"," ")] = line[0].downcase.sub("_"," ")}
end
fuzzy_match_table = FuzzyMatch.new(tpl_lookup.keys)

# Search over all unique names
def search_wrapper(species, tpl_lookup, fuzzy_match_table)
  if tpl_lookup.keys.include? species
    result = species
    distance = -1
  else
    result = fuzzy_match_table.find species
    distance = result.pair_distance_similar species
  end
  return [result, distance]
end
gbif_names = File.readlines(options[:unique]).map {|x| x.chomp.downcase.split(" ")[0..1].join(" ")}
output = Parallel.map(gbif_names, :in_processes => options[:threads], progress: "Fuzzy matching") {|x| search_wrapper(x, tpl_lookup, fuzzy_match_table)}
gbif_lookup = gbif_names.zip(output).to_h

# Write out this lookup file
CSV.open(options[:tpl_cache], "w") do |file|
  file << ["gbif_name", "matched_named", "score"]
  gbif_lookup.each do |spp, matched|
    file << [spp, matched[0], matched[1]]
  end
end

# Do work
CSV.open(options[:output], "w") do |file|
  file << ["species", "match_score", "latitude", "longitude", "year", "match_score"]
  CSV.foreach(options[:gbif], col_sep: "\t").with_index do |line, i|
    if i == 0 then next end
    species = line[0].downcase.split(" ")[0..1].join(" ")
    if gbif_lookup.include? species then puts "#{species}: #{i}" end
    file << [gbif_lookup[species], line[1..line.length]].flatten
  end
end
