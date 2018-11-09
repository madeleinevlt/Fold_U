#! /usr/bin/env Rscript

### FUNCTIONS ###
#################

generate_sum_benchmark_list <- function(benchmark) {
  sum = 0
  sum_benchmark = c()
  family_benchmark = c()
  for (b in benchmark) {
    # A benchmark is found
    if (b != 0){
      sum = sum + 1
      # A family or superfamily benchmark is found
      if (b == 2 | b==3) {
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
  type = scores[,c("benchmark", type_scores)]
  type = type[order(type[,type_scores], decreasing = T),]
  bench = generate_sum_benchmark_list(type$benchmark)
  type = cbind(type, bench)
  return(type)
}

generate_alea_benchmark <- function(N, total_bench) {
  aleas = cbind()
  n = 1
  for (i in 1:n) {
    alea = rep(0, N)
    alea[sample(1:N, total_bench, replace=F)] = 1
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
  return(cbind(scores[,"benchmark"], sum_benchmark))
}


### MAIN ###
############

# Number of templates
N = 412

# Arguments managment
args = commandArgs(trailingOnly=TRUE)
input = paste(args[1], "/scores.csv", sep="")
output = paste(args[1], "/benchmark_rank.png", sep="")

# Read the scores.out file
scores = read.csv(input, header=T)

# Number of benchmark cumulated
ali = create_benchmarch_df(scores, "alignment")
thr = create_benchmarch_df(scores, "threading")
blosum = create_benchmarch_df(scores, "blosum")
scores_sum = create_benchmarch_df(scores, "sum.scores")

# Generate n benchmarks randomly along the ranking
alea = generate_alea_benchmark(N, ali$sum_benchmark[N])

# Create a plot with the rank in x axis and the number of benchmarks in y axis.
# A line representes the benchmarks (family, superfamily and fold) for the different scores generated
# A cross representes the family and superfamily benchmarks

colors = c("black", "red", "orange", "green", "purple")
type_scores = list(alea, ali, thr, blosum, scores_sum)

png(filename=output, height=700, width=700)
  par(mgp=c(2.5,1,0))

  # alignment score
  plot(alea[,"sum_benchmark"], col=colors[1],type='s',
       lwd=3, cex.lab=2.2, cex.axis=2, xlab="Rank", ylab="Benchmark")
  for (i in 2:length(type_scores)) {
    lines(type_scores[[i]]$sum_benchmark, col=colors[i], type="s", lwd=3)
    points(type_scores[[i]]$family_benchmark, pch=4, col=colors[i], cex=2.5, lwd=2)
  }

  legend("bottomright", cex=1.2, lwd=3, col=colors, pch=4,
        c("random", "alignment score", "threading score",
          "blosum score", "alignment + threading + blosum scores"))
dev.off()
