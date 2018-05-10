library(igraph)
library(GroupAssignment)

set.seed(712)

x <- sample_smallworld(1, 20, 4, .1)
d <- distances(x) + 1
d <- max(d) - d + 1

ans <- optimalAssignment(d, leaders = 3, minGroupSize = 2, maxGroupSize = 10)

