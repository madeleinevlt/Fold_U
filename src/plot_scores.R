setwd("~/MeetU/Fold_U/res/His_biosynth")


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

N = 412

scores = read.table("scores.out", header=T)
scores = scores[-1,]

# Normalization
ali_scores = normalization(scores$alignment)
thr_scores = normalization(scores$threading)
blossum_scores  = normalization(scores$blossum)
thr_blossum_scores =  normalization(thr_scores + blossum_scores)
ali_thr_blossum_scores =  normalization(ali_scores + thr_scores + blossum_scores)


ali = cbind(scores[,c(1,2)], ali_scores)
ali = ali[order(ali$ali_scores, decreasing = T),]
bench = generate_sum_benchmark_list(ali$benchmark)
ali = cbind(ali, bench)

thr = cbind(scores[,c(1,2)], thr_scores)
thr = thr[order(thr$thr_scores, decreasing = T),]
bench = generate_sum_benchmark_list(thr$benchmark)
thr = cbind(thr, bench)

blossum = cbind(scores[,c(1,2)], blossum_scores)
blossum = blossum[order(blossum$blossum_scores, decreasing = T),]
bench = generate_sum_benchmark_list(blossum$benchmark)
blossum = cbind(blossum, bench)

thr_blossum = cbind(scores[,c(1,2)], thr_blossum_scores)
thr_blossum = thr_blossum[order(thr_blossum$thr_blossum_scores, decreasing = T),]
bench = generate_sum_benchmark_list(thr_blossum$benchmark)
thr_blossum = cbind(thr_blossum, bench)

ali_thr_blossum = cbind(scores[,c(1,2)], ali_thr_blossum_scores)
ali_thr_blossum = ali_thr_blossum[order(ali_thr_blossum$ali_thr_blossum_scores, decreasing = T),]
bench = generate_sum_benchmark_list(ali_thr_blossum$benchmark)
ali_thr_blossum = cbind(ali_thr_blossum, bench)


png(filename="benchmark_rank.png", height=700, width=700)
par(mgp=c(2.5,1,0))
# alignment score
plot(ali$sum_benchmark, col="green",type='s', lwd=3, cex.lab=2.2, cex.axis=2, xlab="Rank", ylab="Benchmark")
points(ali$family_benchmark, pch=4, col="green", cex=2.5, lwd=2)

# threading score
lines(thr$sum_benchmark, col="red", type="s", lwd=3)
points(thr$family_benchmark, pch=4, col="red", cex=2.5, lwd=2)

# blossum score
lines(blossum$sum_benchmark, col="blue", type="s", lwd=3)
points(blossum$family_benchmark, pch=4, col="blue", cex=2.5, lwd=2)

# threading + blossum scores
lines(thr_blossum$sum_benchmark, col="orange", type="s", lwd=3)
points(thr_blossum$family_benchmark, pch=4, col="orange", cex=2.5, lwd=2)

# alignment + threading + blossum scores
lines(ali_thr_blossum$sum_benchmark, col="green", type="s", lwd=3)
points(ali_thr_blossum$family_benchmark, pch=4, col="green", cex=2.5, lwd=2)


abline(b=ali$sum_benchmark[N]/N, a=0, lwd=3)
legend("bottomright", cex=0.8, c("alignment score", "threading score", "blossum score", "threading + blossum scores", "alignment + threading + blossum scores"), lwd=3, col=c("green", "red", "blue", "orange", "yellow"), pch=4)

dev.off()

# png(filename="legend.png", height=700, width=700)
# plot(c(0,10),c(0,10))
# legend(c(5,5), cex=2, c("alignment score", "threading score", "threading + blossum scores"), lwd=3, col=c("green", "red", "blue"), pch=4)
# dev.off()
