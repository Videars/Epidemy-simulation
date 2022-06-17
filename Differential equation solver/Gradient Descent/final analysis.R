library(igraph)
library(ggplot2)
library(deSolve)
library(reshape2)
library(R.matlab)
library(Matrix)
library(poweRlaw)

sol_hom <- unlist(readMat('sol_hom.mat'))
sol_sm <- unlist(readMat('sol_sm.mat'))
sol_sf <- unlist(readMat('sol_sf.mat'))

sol_hom_mean = mean(sol_hom)
sol_sm_mean = mean(sol_sm)
sol_sf_mean = mean(sol_sf)

a_hom <- matrix(unlist(readMat('A_hom.mat')), 100,100)
a_sm <- matrix(unlist(readMat('A_sm.mat')), 100,100)
a_sf <- matrix(unlist(readMat('A_sf.mat')), 100,100)

hom_graph <- graph_from_adjacency_matrix(a_hom)
sm_graph <- graph_from_adjacency_matrix(a_sm)
sf_graph <- graph_from_adjacency_matrix(a_sf)

hom_degree <- degree(hom_graph )
sm_degree <- degree(sm_graph)
sf_degree <- degree(sf_graph)

hom_borders <- c(min(sol_hom), max(sol_hom))
sm_borders <- c(min(sol_sm), max(sol_sm))
sf_borders <- c(min(sol_sf), max(sol_sf))

plot(hom_degree, sol_hom, xlab="grado", ylab="punto di stabilità", main="punti di stabilità su grado per Omogeneo")
abline(h=seq(min(sol_hom),max(sol_hom),0.01), col="gray", lty="dotted")
abline(h=sol_hom_mean, col='red')

plot(sm_degree, sol_sm, xlab="grado", ylab="punto di stabilità", main="punti di stabilità su grado per Small World")
abline(h=seq(min(sol_sm),max(sol_sm),0.01), col="gray", lty="dotted")
abline(h=sol_sm_mean, col='red')

plot(sf_degree, sol_sf, xlab="grado", ylab="punto di stabilità", main="punti di stabilità su grado per Scale Free")
abline(h=seq(min(sol_sf),max(sol_sf),0.01), col="gray", lty="dotted")
abline(h=sol_sf_mean, col='red')