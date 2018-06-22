#!/usr/bin/python
#Parsing NCBI Taxonomy to pull out all plant names
#Will Pearse - 2013-09-11

#Set the node we're interested in - green plants
head_node = 33090
nodes = set([head_node])
species = set([-1])

#Process the list given to it
# - I've checked, and deleting like this is OK. Which is strange, because I thought it wasn't...
def process_list(node_file, nodes, max_number):
    for i,node_line in enumerate(node_file):
        if i == max_number:
            return (node_file, nodes, species)
        split = node_line.split("\t|\t")
        if int(split[1]) in nodes:
            nodes.add(int(split[0]))
            if split[2] == "species":
                species.add(int(split[0]))
            del node_file[i]
    return (node_file, nodes, species)

#Loop over the subsets of the data
node_file = open("nodes.dmp").readlines()
for i in range(1,101):
    print i, int(len(node_file) * (i/100.0)), len(nodes), len(species)
    i = int(len(node_file) * (i/100.0))
    node_file, nodes, species = process_list(node_file, nodes, i)

#Now iterate until we have everything
last_val = len(nodes)
i = 1
while last_val != len(nodes):
    node_file, nodes, species = process_list(node_file, nodes, len(node_files)+100)
    print len(nodes), len(species), i
    i += 1

#Grab TaxonIDs and names
# - ignoring synonyms...
names = []
with open("names.dmp") as name_file:
    for name_line in name_file:
        split = name_line.split("\t|\t")
        split[3] = split[3].strip().replace("\t|", "")
        if int(split[0]) in species and split[3] == "scientific name":
            names.append(split[1])

with open("genbank_raw_names.txt", "w") as names_output:
    for each in names:
        names_output.write(str(each) + "\n")
