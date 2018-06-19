require 'rake'
require 'set'
require 'unidecoder'

task :default => [:build]

desc "Run analyses"
task :build => [:download, :harmonise]

################################
# Tidying up ###################
################################
require 'rake/clean'
CLOBBER.include("raw_data/*")
CLEAN.include("clean_data/*")

################################
# Download raw data ############
################################
# Wrapper
desc "Download all raw data"
task :download => [:dwn_setup, :dwn_genbank, :dwn_try, :dwn_gbif]

# Setup
desc "Download setup"
task :dwn_setup do
  unless File.directory?("raw_data") then Dir.mkdir("raw_data") end
end

# GenBank
def dwn_genbank
  Dir.chdir("raw_data") do 
    `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`
    `tar -xf taxdump.tar.gz nodes.dmp`
    `tar -xf taxdump.tar.gz names.dmp`
    File.delete("taxdump.tar.gz")
  end
end
file 'raw_data/names.dmp' do
  dwn_genbank
end
file 'raw_data/nodes.dmp' do
  dwn_genbank
end
desc "Download GenBank data"
task :dwn_genbank => ["raw_data/names.dmp", "raw_data/nodes.dmp"]

# TRY
file 'raw_data/TryAccSpecies.txt' do
  Dir.chdir("raw_data") do 
    `wget https://www.try-db.org/dnld/TryAccSpecies.txt --no-check-certificate`
  end
end
desc "Download TRY data"
task :dwn_try => "raw_data/TryAccSpecies.txt"

# GBIF
file "raw_data/gbif_cut.txt", [:location] do |task, args|
  Dir.chdir("raw_data") do
    if args[:location].nil?
      location = "http://api.gbif.org/v1/occurrence/download/request/0004202-180508205500799.zip"
      puts "Using default GBIF data url: #{args[:location]}"
      puts "To change, run with something like 'rake raw_data/gbif_cut.txt[http://api.gbif.org/url/to/gbif/00001234.zip]"
    else
      location = args[:location]
    end
    file_stem = location[/[0-9\\-]{10,}/]
    `wget #{location}`
    `unzip #{file_stem}.zip`
    `cut -f13,17,18,28 #{file_stem}.csv > gbif_cut.txt`
    File.delete("#{file_stem}.zip")
    File.delete("#{file_stem}.csv")
  end
end
desc "Download GBIF data"
task :dwn_gbif => "raw_data/gbif_cut.txt"

# GIS layers
def dwn_gis
  Dir.chdir("raw_data") do 
    `wget http://legacy.jefferson.kctcs.edu/techcenter/gis%20data/World/Zip/Continents.zip`
    `unzip Continents.zip`
    File.delete("Continents.zip")
  end
end
file 'raw_data/continent.avl' do dwn_gis end
file 'raw_data/continent.dbf' do dwn_gis end
file 'raw_data/continent.htm' do dwn_gis end
file 'raw_data/continent.prj' do dwn_gis end
file 'raw_data/continent.sbn' do dwn_gis end
file 'raw_data/continent.shp' do dwn_gis end
file 'raw_data/continent.shp.xml' do dwn_gis end
file 'raw_data/continent.shx' do dwn_gis end
file 'raw_data/Continents.lyr' do dwn_gis end
file 'raw_data/Continents.lyr.xml' do dwn_gis end
desc "Download GIS data"
task :dwn_gis => ["raw_data/continent.avl","raw_data/continent.dbf","raw_data/continent.htm","raw_data/continent.prj","raw_data/continent.sbn","raw_data/continent.shp","raw_data/continent.shp.xml","raw_data/continent.shx","raw_data/Continents.lyr","raw_data/Continents.lyr.xml"]

################################
# Taxonomy cleaning ############
################################
# Wrapper
desc "Harmonise all raw data"
task :harmonise => [:hrm_setup, :hrm_genbank, :hrm_try, :hrm_gbif]  

# Setup
desc "Clean GBIF data"
task :hrm_setup do
  unless File.directory?("clean_data") then Dir.mkdir("clean_data") end
end

# GBIF
file 'clean_data/gbif_spp_clean.txt' do
  puts "GBIF doesn't have authors etc. cut so this isn't perfect"
  Dir.chdir("raw_data") do
    all_spp = Set.new
    File.open("gbif_cut.txt", "r") do |file|
      file.each_with_index do |line, i|
        if i == 0 then next end
        line = line.to_ascii
        line = line.split("\t")
        species = line[0]
        species.downcase!
        all_spp.add species
      end
    end
    File.open("../clean_data/gbif_spp_clean.txt", "w") {|file| all_spp.each {|x| file << "#{x}\n"}}
  end
end
desc "Clean GBIF data"
task :hrm_gbif => "clean_data/gbif_spp_clean.txt"

# TRY
file 'clean_data/try_spp_clean.txt' do
  Dir.chdir("raw_data") do
    puts "TRY hybrids are cut at |Genus x| genus, so this isn't perfect"
    all_spp = Set.new
    File.open("TryAccSpecies.txt", "r") do |file|
      file.each_with_index do |line, i|
        if i == 0 then next end
        # The blasted hybrid sign
        begin
          line = line.to_ascii
          line = line.split "\t"
        rescue ArgumentError
          line = line.encode("UTF-8", "Windows-1252")
          line = line.sub("\u00D7", "x ")
          line.sub("x  ", "x ")
          line = line.to_ascii
          line = line.split "\t"          
        end
        species = line[1].gsub('"', "")
        species.downcase!
        species = species.split(" ")[0..1].join(" ")
        all_spp.add species
      end
    end
    File.open("../clean_data/try_spp_clean.txt", "w") {|file| all_spp.each {|x| file << "#{x}\n"}}
  end
end
desc "Clean TRY data"
task :hrm_try => "clean_data/try_spp_clean.txt"

# GenBank
file 'clean_data/genbank_spp_clean.txt' do
  Dir.chdir("raw_data") do
    `../src/parse_genbank.py`
    all_spp = Set.new
    File.open("genbank_raw_names.txt", "r") do |file|
      file.each_with_index do |line, i|
        species = line.downcase
        species = species.to_ascii
        all_spp.add species
      end
    end
    File.open("../clean_data/genbank_spp_clean.txt", "w") {|file| all_spp.each {|x| file << "#{x}\n"}}
  end
end
desc "Clean GenBank data"
task :hrm_genbank => "clean_data/genbank_spp_clean.txt"
