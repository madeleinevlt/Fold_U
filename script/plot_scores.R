#! /usr/bin/env Rscript

normalization <- function(scores) {
  return((scores - min(scores)) / (max(scores) - min(scores)))
}

generate_sum_benchmark_list <- function(benchmark) {
  sum = 0
  sum_benchmark = c()
  family_benchmark = c()
  for (b in benchmark) {
    if (b != "."){
      sum = sum + 1
      if (b == "Family" | b=="Superfamily") {
        family_benchmark = c(family_benchmark, sum)
      } else {
        family_benchmark = c(family_benchmark, -10)
      }
    } else {
      family_benchmark = c(family_benchmark, -10)
    }
    sum_benchmark = c(sum_benchmark, sum)
  }
  return(cbind(sum_benchmark, family_benchmark))
}

create_benchmarch_df <- function(scores, type_scores) {
  type = cbind(scores[,c(1,2)], type_scores)
  type = type[order(type$type_scores, decreasing = T),]
  bench = generate_sum_benchmark_list(type$benchmark)
  type = cbind(type, bench)
  return(type)
}

N = 412

# Arguments managment
args = commandArgs(trailingOnly=TRUE)
input = paste(args[1], "/scores.out", sep="")
output = paste(args[1], "/benchmark_rank.png", sep="")

# Reading of the scores.out file
scores = read.table(input, header=T)
scores = scores[-1,]

# Normalization of the scores using min-max scaling method
ali_scores = normalization(scores$alignment)
thr_scores = normalization(scores$threading)
blosum_scores  = normalization(scores$blosum)
ali_thr_scores = normalization(ali_scores + thr_scores)
thr_blosum_scores =  normalization(thr_scores + blosum_scores)
ali_thr_blosum_scores =  normalization(ali_scores + thr_scores + blosum_scores)

# Number of benchmark cumulated
ali = create_benchmarch_df(scores, ali_scores)
thr = create_benchmarch_df(scores, thr_scores)
blosum = create_benchmarch_df(scores, blosum_scores)
ali_thr = create_benchmarch_df(scores, ali_thr_scores)
thr_blosum = create_benchmarch_df(scores, thr_blosum_scores)
ali_thr_blosum = create_benchmarch_df(scores, ali_thr_blosum_scores)

# Generate n benchmarks randomly along the ranking
aleas = cbind()
n = 1
for (i in 1:n) {
  alea = rep(0, N)
  alea[sample(1:N, ali$sum_benchmark[N], replace=F)] = 1
  sum = 0
  sum_benchmark = c()
  for (b in alea) {
    if (b == 1) {
      sum = sum + 1
    }
    sum_benchmark = c(sum_benchmark, sum)
  }
 aleas = cbind(aleas,sum_benchmark)
}
sum_benchmark = rowSums(aleas)/n
alea = cbind(scores[,c(1,2)], sum_benchmark)

# Create a plot with the rank in x axis and the number of benchmarks in y axis.
# A line representes the benchmarks (family, superfamily and fold) for the different scores generated
# A cross representes the family and superfamily benchmarks

colors = c("black", "red", "orange", "green", "purple")
type_scores = list(alea, ali, thr, blosum, ali_thr_blosum)

png(filename=output, height=700, width=700)
  par(mgp=c(2.5,1,0))

  # alignment score
  plot(type_scores[[1]]$sum_benchmark, col=colors[1],type='s',
       lwd=3, cex.lab=2.2, cex.axis=2, xlab="Rank", ylab="Benchmark")
  for (i in 2:length(type_scores)) {
    lines(type_scores[[i]]$sum_benchmark, col=colors[i], type="s", lwd=3)
    points(type_scores[[i]]$family_benchmark, pch=4, col=colors[i], cex=2.5, lwd=2)
  }

  legend("bottomright", cex=1.2, lwd=3, col=colors, pch=4,
        c("random", "alignment score", "threading score",
          "blosum score", "alignment + threading + blosum scores"))
dev.off()
