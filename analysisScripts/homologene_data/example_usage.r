source('homologene_data/homoLookup.r')
homoDat <- homoLookup(9606,10116) # this would be for human and rat

# this returns a lookup table with the homogene, the entrez id in species one and entrez id in species two, thereby allowing you to lookup what are the homologs of a particular gene.
