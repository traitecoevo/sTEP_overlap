require 'rake'
require 'set'
require 'unidecoder'
require 'csv'
require 'os'

# Wrapper for line splitting with hybrid problem
def line_split(line, split="\t")
  begin
    line = line.to_ascii
    line = line.split(split)
  rescue ArgumentError
    line = line.encode("UTF-8", "Windows-1252")
    line.sub!("\u00D7", "x ")
    line.sub!("x  ", "x ")
    line = line.to_ascii
    line = line.split(split)
  end
  return line
end

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
task :download => [:dwn_setup, :dwn_tpl, :dwn_try, :dwn_genbank, :dwn_gbif]

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
file 'raw_data/names.dmp' do dwn_genbank end
file 'raw_data/nodes.dmp' do dwn_genbank end
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
    `7z e #{file_stem}.zip`
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
    `wget http://faculty.baruch.cuny.edu/geoportal/data/esri/world/continent.zip`
    `7z e continent.zip`
    File.delete("continent.zip")
  end
end
file 'raw_data/continent.dbf' do dwn_gis end
file 'raw_data/continent.prj' do dwn_gis end
file 'raw_data/continent.sbn' do dwn_gis end
file 'raw_data/continent.sbx' do dwn_gis end
file 'raw_data/continent.shp' do dwn_gis end
file 'raw_data/continent.shp.xml' do dwn_gis end
file 'raw_data/continent.shx' do dwn_gis end
desc "Download GIS data"
task :dwn_gis => ["raw_data/continent.dbf","raw_data/continent.prj","raw_data/continent.sbn","raw_data/continent.sbx","raw_data/continent.shp","raw_data/continent.shp.xml","raw_data/continent.shx"]

# TPL species lists
file 'raw_data/tpl_names.txt' do
  Dir.chdir("raw_data") do 
    `wget https://raw.githubusercontent.com/schwilklab/taxon-name-utils/master/data/theplantlist1.1/names_unique.csv`
    File.rename("names_unique.csv", "tpl_names.txt")
    if OS.linux?
      `tr '\r' '\n' < tpl_names.txt > tmp.txt`
      `mv tmp.txt tpl_names.txt`
    end
  end
end
desc "Download The Plant List accepted names"
task :dwn_tpl => "raw_data/tpl_names.txt"


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
def hrm_gbif
  Dir.chdir("raw_data") do
    # Build TPL set
    tpl = Set.new
    File.open("tpl_names.txt") do |file|
      file.each_with_index do |line, i|
        if i == 0 then next end
        sp = line_split(line, ",")[5]
        sp = sp.sub("_", " ").downcase.chomp
        tpl.add(sp)
      end
    end
    # Parse GBIF data
    all_spp = Set.new
    File.open("gbif_cut.txt", "r") do |file|
      File.open("../clean_data/gbif_tpl_locations.txt", "w") do |spp_tpl|
        file.each_with_index do |line, i|
          # First handle file loading (can't use CSV because of encoding)
          if i == 0 then next end
          if i < 10000 then next end
          line = line_split(line)
          # Make species set
          species = line[0]
          species.downcase!
          all_spp.add species
          # Make TPL-checked lats and longs (as we go along)
          species = species.split(" ")[0..1].join " "
          if tpl.include? species
            line[0] = species
            spp_tpl << line.join("\t")
          end
        end
      end
    end
    # Write out all GBIF species (at the end, from the set)
    File.open("../clean_data/gbif_spp.txt", "w") {|file| all_spp.each {|x| file << "#{x}\n"}}
  end
end
file 'clean_data/gbif_spp.txt' do hrm_gbif end
file 'clean_data/gbif_tpl_locations.txt' do hrm_gbif end
desc "Clean GBIF data"
task :hrm_gbif => ["clean_data/gbif_spp.txt", "clean_data/gbif_tpl_locations.txt"]

# TRY
file 'clean_data/try_spp_clean.txt' do
  Dir.chdir("raw_data") do
    puts "TRY hybrids are cut at |Genus x| genus, so this isn't perfect"
    all_spp = Set.new
    File.open("TryAccSpecies.txt", "r") do |file|
      file.each_with_index do |line, i|
        if i == 0 then next end
        line = line_split(line)
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

################################
# Analysis #####################
################################
# Wrapper
desc "Run analyses"
task :analysis => [:analysis_setup, 'figures/multi_gam.png', 'figures/order_phy.pdf', 'figures/multi_gam.png', 'tables/two_and_three_way_comparisons.tex', 'tables/all_families_ranking.csv', 'tables/summary_of_endemic_analysis.csv', :analysis_folder]

# Setup
desc "Install R packages"
task :analysis_setup do
  `Rscript R/install.R`
  unless File.directory?("output") then Dir.mkdir("output") end
  unless File.directory?("overlap_data") then Dir.mkdir("overlap_data") end
  unless File.directory?("tables") then Dir.mkdir("tables") end
  unless File.directory?("figures") then Dir.mkdir("figures") end
end

# Getting Things Done
file 'figures/multi_gam.png' do `Rscript -e "source('R/install.R)'; run_gam_df()"` end
file 'figures/order_phy.pdf' do `Rscript -e "source('R/install.R)'; makephylo()"` end
file 'figures/multi_gam.png' do `Rscript -e "source('R/install.R'); multi_gam()"` end
file 'tables/two_and_three_way_comparisons.tex' do `Rscript -e "source(R/install.R); write_overlap_table(overlap_data)"` end
file 'tables/all_families_ranking.csv' do `Rscript -e "source('R/install.R'); do_big_list_family_anlysis()")` end
file 'tables/summary_of_endemic_analysis.csv' do `Rscript -e "source('R/install.R'); do.endemic.analysis()"` end
task :analysis_folder do
  `Rscript -e "source('R/install.R'); do_overlap_analysis()"`
end
