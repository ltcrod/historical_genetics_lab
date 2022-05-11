#! /usr/bin/env Rscript
library(reshape)

args <- commandArgs(trailing = TRUE)
fn <- args[1]
fn_in = paste(fn, "_relatedness_per_individual.txt", sep="")
fn_out = paste(fn, "_relatedness_per_individual.wide.txt", sep="")

long=read.csv(fn_in,header=T, sep='\t')
wid=reshape(long, idvar=c("id"), timevar = "relatedness", direction="wide")
write.table(wid, fn_out, sep='\t', quote=F, row.names=F)

