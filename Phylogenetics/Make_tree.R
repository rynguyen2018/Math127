library(ggtree)
nwk<- "/Users/ryannguyen/Desktop/Math_127/Phylogenetics/tree2.nwk"
tree<- read.tree(nwk)

ggtree(tree, color="firebrick", lbranch.length='none') + geom_tiplab() +geom_tippoint()

ggsave("plswork.png")

#0.1
#0.1
#0.1
#0.3
#0.3
#0.3
#0.3
#0.5
#0.2
#0.2
#0.2