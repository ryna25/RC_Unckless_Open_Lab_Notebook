#Set Working Directory
setwd("/Users/anjaligupta/Google Drive/My Drive/UncklessLab/Personnel/Anjali/Computational_Pipelines")

#Load required packages
library(ape)
library(phangorn)
library(spatstat)
library(spatstat.geom)
library(phytools)

##Import output tree file

#First I'm doing this for the X

x <- phytools::read.newick("OutputTrees_X_100kb.trees")

#Plot overlaid gene trees

##root each one using pseudoobscura
for(i in 1:30){
  x[[i]]<-root(x[[i]],"Dpse.bam")
}

#
##make each one ultrametric

for(i in 1:30){
  x[[i]]<-chronos(x[[i]])
}

#plot overlaid
#densiTree(x, type="cladogram",nodes="intermediate", alpha = .1, consensus = x[[1]]$tip.label[c(11,6:10,12:15,1:5)])
densiTree(x, type="cladogram",nodes="intermediate", alpha = .05, consensus = x[[1]]$tip.label[c(1:5,6:9,13:14,12,10,11,15)])

#manually move the branch leading to the outgroup behind the tree for easier viewing
for(i in 1:30){
  x[[i]]$edge.length[x[[i]]$edge.length==1]<- -10
}

#replot
densiTree(x, type="cladogram",nodes="intermediate", alpha = .2, consensus = x[[1]]$tip.label[c(1:5,6:9,13:14,12,10,11,15)])

#how often is affinis monophyletic
x.mono<-c()
for(i in 1:30){
  x.mono[i]<-(is.monophyletic(x[[i]],x[[1]]$tip.label[c(1:3)]))
}
table(x.mono)






#Now, I'm doing this for the autosomes
Autosome <- phytools::read.newick("OutputTrees_Autosome_100kb.trees")

#Plot overlaid gene trees

##root each one using pseudoobscura
for(i in 1:60){
  Autosome[[i]]<-root(Autosome[[i]],"Dpse.bam")
}

#
##make each one ultrametric

for(i in 1:60){
  Autosome[[i]]<-chronos(Autosome[[i]])
}

#plot overlaid
#densiTree(Autosome, type="cladogram",nodes="intermediate", alpha = .1, consensus = x[[1]]$tip.label[c(11,6:10,12:15,1:5)])
densiTree(Autosome, type="cladogram",nodes="intermediate", alpha = .05, consensus = x[[1]]$tip.label[c(1:5,6:9,13:14,12,10,11,15)])

#manually move the branch leading to the outgroup behind the tree for easier viewing
for(i in 1:60){
  Autosome[[i]]$edge.length[Autosome[[i]]$edge.length==1]<- -10
}

#replot
densiTree(Autosome, type="cladogram",nodes="intermediate", alpha = .2, consensus = x[[1]]$tip.label[c(1:5,6:9,13:14,12,10,11,15)])


#how often is affinis monophyletic
Autosome.mono<-c()
for(i in 1:60){
  Autosome.mono[i]<-(is.monophyletic(Autosome[[i]],Autosome[[1]]$tip.label[c(1:3)]))
}
table(Autosome.mono)
