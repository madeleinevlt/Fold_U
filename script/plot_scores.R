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
intput = paste(args[1], "/scores.out", sep="")
output = paste(args[1], "/benchmark_rank.png", sep="")

# Reading of the scores.out file
scores = read.table(intput, header=T)
scores = scores[-1,]

# Normalization of the scores using min-max scaling method
ali_scores = normalization(scores$alignment)
thr_scores = normalization(scores$threading)
blossum_scores  = normalization(scores$blossum)
ali_thr_scores = normalization(ali_scores + thr_scores)
thr_blossum_scores =  normalization(thr_scores + blossum_scores)
ali_thr_blossum_scores =  normalization(ali_scores + thr_scores + blossum_scores)

# Number of benchmark cumulated
ali = create_benchmarch_df(scores, ali_scores)
thr = create_benchmarch_df(scores, thr_scores)
blossum = create_benchmarch_df(scores, blossum_scores)
ali_thr = create_benchmarch_df(scores, ali_thr_scores)
thr_blossum = create_benchmarch_df(scores, thr_blossum_scores)
ali_thr_blossum = create_benchmarch_df(scores, ali_thr_blossum_scores)


# Create a plot with the rank in x axis and the number of benchmarks in y axis.
# A line representes the benchmarks (family, superfamily and fold) for the different scores generated
# A cross representes the family and superfamily benchmarks

colors = c("red", "orange", "green", "purple")
type_scores = list(ali, thr, blossum, ali_thr_blossum)

png(filename=output, height=700, width=700)
  par(mgp=c(2.5,1,0))

  # alignment score
  plot(type_scores[[1]]$sum_benchmark, col=colors[1],type='s', lwd=3, cex.lab=2.2, cex.axis=2, xlab="Rank", ylab="Benchmark")
  points(type_scores[[1]]$family_benchmark, pch=4, col=colors[1], cex=2.5, lwd=2)
  abline(b=ali$sum_benchmark[N]/N, a=0, lwd=3)

  for (i in 2:length(type_scores)) {
    lines(type_scores[[i]]$sum_benchmark, col=colors[i], type="s", lwd=3)
    points(type_scores[[i]]$family_benchmark, pch=4, col=colors[i], cex=2.5, lwd=2)
  }

  legend("bottomright", cex=1.2, lwd=3, col=colors, pch=4,
        c("alignment score", "threading score", "blossum score",
          "alignment + threading + blossum scores"))
dev.off()
