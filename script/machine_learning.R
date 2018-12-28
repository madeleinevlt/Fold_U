# Install packages if not already installed and load them
list.of.packages <- c("boot", "dplyr", "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(boot)
library(dplyr)
library(readr)

# Reduce verbose of readr package when readin csv files
options(readr.num_columns = 0)

prots.names = list.dirs(path = "./results", full.names = FALSE, recursive = FALSE)
path.to.csv.scores = paste0("./results/", prots.names, "/scores.csv")
path.to.csv.scores = path.to.csv.scores[which(path.to.csv.scores != "./results/plots/scores.csv")]

# Read the result scores into a dataframe
benchmark <- path.to.csv.scores %>% 
             lapply(read_csv) %>% 
             bind_rows
df.benchmark = as.data.frame(benchmark)

# Visualize the distribution of scores values with boxplots
png(file="./results/plots/scores_distribution.png", units="px", width=1600, height=1600, res=300)
par(mar=c(10, 4, 4, 2) + 0.1)
a = boxplot(df.benchmark[,-c(1, 2)], main="Distribution of values for each scores", xlab="", ylab="Values", las=2, col=rainbow(7), xaxt="n")
mtext(text="Scores", side=1, line=7)
labs = colnames(df.benchmark[-c(1,2)])
text(seq_along(labs), par("usr")[3], labels = labs, srt = 35, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
dev.off()

# Transform the benchmarks (".", "Family", "Superfamily", "Fold") into binary numeric values:
# "Family", "Superfamily", "Fold" -> 1
#                             "." -> 0
new.Y = rep(0, dim(df.benchmark)[1])
new.Y[which(df.benchmark[, (which(colnames(df.benchmark) == "benchmark"))] != ".")] = 1

# Build the data.frame with benchmarks as factors (0, 1)
df.x = df.benchmark[, -(which(colnames(df.benchmark) == "benchmark"))]
df.XY = data.frame(df.x, as.factor(new.Y))

# Transform the NAs -> 0
df.benchmark[is.na(df.benchmark)] = 0

# Set column names of the data.frame
colnames(df.XY) = c(colnames(df.x),"benchmark")

# Creation of the training and validation samples: 2/3 (train) - 1/3 (validate)
v.ind.train = sample(1:nrow(df.XY), size = 2/3 * nrow(df.XY), replace=FALSE)
v.ind.test = (1:nrow(df.XY))[-v.ind.train]
train.mat = df.XY[v.ind.train,]
test.mat = df.XY[v.ind.test,]
# Clean test and train matrices
train.mat = train.mat[,-c(1, 8)]
test.mat = test.mat[,-c(1, 8)]


vIndApp=sample(1:nrow(df.XY),size=2/3*nrow(df.XY),replace=FALSE)
vIndTest=(1:nrow(df.XY))[-vIndApp]
    
matApp=df.XY[vIndApp,]
matTest=df.XY[vIndTest,]


# Launch Generalized Linear Model
# family = binomial: because benchmark is binary 0 or 1
fit = glm(benchmark ~ ., data=train.mat, family="binomial")
summary(fit)
# Retrieve the new scores coefficients
coefficients = exp(coef(fit))
cat("Coefficients (weights) :", coefficients,"\n")

# Write new scores.csv adding a column for the weighted score
for (scores.file in path.to.csv.scores){
    tmp <- read.csv(scores.file)
    tmp <- as.data.frame(tmp)
    # Set the new weighted summed scores based on the calculated coefficients as:
    weighted_combined_scores =  coefficients[1] +
                                coefficients[2] * tmp$alignment +
                                coefficients[3] * tmp$threading +
                                coefficients[4] * tmp$modeller +
                                coefficients[5] * tmp$secondary_structure +
                                coefficients[6] * tmp$solvent_access +
                                coefficients[7] * tmp$co_evolution
    # Normalisation with min-max Z-Score: 0 <= score <= 1 
    weighted_combined_scores = (weighted_combined_scores - min(weighted_combined_scores)) / (max(weighted_combined_scores) - min(weighted_combined_scores))
    tmp <- cbind(tmp, weighted_combined_scores)
    # Write the CSV file
    write.csv(tmp, scores.file, row.names = F, quote = F)
}


##### Robust Z-Score for normalisation is not giving better results:
    # m <- median(df.benchmark$sum_scores)
    # s <- mad(df.benchmark$sum_scores)
    # df.benchmark$sum_scores <- abs(df.benchmark$sum_scores - m) / (s)
#####


# Calculate the performances of the model
yobsApp=train.mat$benchmark
ppredApp=predict(fit, newdata=train.mat, type = "response")

yobsTest=test.mat$benchmark
ppredTest=predict(fit, newdata=test.mat, type = "response")

yPredApp = round(ppredApp)
yPredTest = round(ppredTest)

tableApp = table(yobsApp,yPredApp)
tableTest = table(yobsTest, yPredTest)

# TRAINING
# vp/vp+fn
sensiA = tableApp[4] / (tableApp[4] + tableApp[2])

# vn/vn+fp
speciA = tableApp[1] / (tableApp[1] + tableApp[3])

# vp+vn/tot
tauxA = (tableApp[4] + tableApp[1]) / (tableApp[1] + tableApp[2] + tableApp[3] + tableApp[4])

# TEST
# vp/vp+fn
sensiP = tableTest[4] / (tableTest[4] + tableTest[2])

# vn/vn+fp
speciP = tableTest[1] / (tableTest[1] + tableTest[3])

# vp+vn/tot
tauxP = (tableTest[4] + tableTest[1]) / (tableTest[1] + tableTest[2] + tableTest[3] + tableTest[4])


# Cross Validation --------------------------------------------------------

### leave-one-out ###
# The first component is the raw cross-validation estimate of prediction error.
# The second component is the adjusted cross-validation estimate.
# The adjustment is designed to compensate for the bias introduced by not using
# leave-one-out cross-validation.
cv.err <- cv.glm(test.mat, fit)$delta
cat("Leave-one-out CV error estimate: ", cv.err[1], "\n")