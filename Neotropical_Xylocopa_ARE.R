################################ Ancestral range estimation ################################
################################       BioGeoBEARS          ################################
############################################################################################


####drop outgroup tips in tree####
library(ape)
tr <- read.nexus("Xylocopa9820168.tre")
tr$tip.label
plot(tr, no.margin=T, cex=0.5)
outgroup <- tr$tip.label[c(58, 4, 98, 14, 9, 10, 3, 89, 1, 94, 95, 13, 18, 58, 59)]
outgroup
tr2 <- drop.tip(tr, outgroup)
plot(tr2, no.margin=T, cex=0.5)
is.binary(tr2)
write.tree(tr2, file="Xylocopa_outgroup_less.tre") #for newick format file
write.nexus(tr2, file="Xylocopa_outgroup_less.tre")
