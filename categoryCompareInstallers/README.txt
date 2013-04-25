These are the source (*.tar.gz) and zip (*.zip) installers for the 
version of categoryCompare used to generate the results
reported in the associated manuscript. 

This version was compiled using R 2.14.1, and should work fine in that
R version. Due to NAMESPACE considerations, and changes to dependent packages,
it may not install properly in newer versions of R. 

Prior to installation, you should first install the major dependencies

biocLite(c("RCytoscape", "Category", "GO.db", "KEGG.db", "colorspace", "graph"))

To redo the analysis reported in the publication you will also need the chip and organism specific packages:

biocLite(c("org.Rn.eg.db", "org.Mm.eg.db", "rat2302.db", "mouse4302.db"))

And some analysis stuff:

biocLite(c("genefilter", "limma", "affy"))

"other_packages" contains two packages that may be helpful when trying to install categoryCompare_0.6.20 in an older version of R: an archived version of DBI that installs in R 2.14.1, and a newer version of RBGL that may be necessary for installation on newer machines, due to changes in the Boost library.
