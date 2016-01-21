# Compiled tree traits
setwd("~/Documents/git/treetraits/")
d <- read.csv("Summer 2015 Tree Traits.csv")

length(unique(d$Species))


dt <- with(d, table(Species, Site))

# Species we have at only one site

atonesite <- apply(dt, 1, function(x) length(x[x==0])==3)

atonesite[atonesite==TRUE]

# At all four sites

atfoursite <- apply(dt, 1, function(x) length(x[x==0])==0)

length(atfoursite[atfoursite==TRUE])


