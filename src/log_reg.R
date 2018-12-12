# Bash

# cd Results
# echo ",benchmark,alignment,threading,modeller,secondary_structure,solvent_access,co_evolution,sum_scores" > tous_les_scores.csv
# for dossier in ./* ; do if [ -d $dossier ]; then tail -q -n +2 $dossier/scores.csv >> tous_les_scores.csv; fi done

# Names of the proteins belonging to the 4 different SCOP classes
# SCOP_A = c("DEP", "DnaJ", "hemery", "TBCA", "histone", "GA")
# SCOP_B = c("igvar-h", "Lum_binding", "Cohesin", "Agglutinin", "SSB")
# SCOP_C = c("His_biosynth", "Rib_hydrolayse", "LRR", "Lipoprotein_4", "ETF_alpha")
# SCOP_D = c("UBQ", "PAC", "svmp", "FAD", "PCNA")

# Read the result scores for every classes of SCOP ids
benchmark = read.table("SCOP_A/scores_scop_class_a.csv", sep=",", header = TRUE)
benchmark = read.table("SCOP_B/scores_scop_class_b.csv", sep=",", header = TRUE)
benchmark = read.table("SCOP_C/scores_scop_class_c.csv", sep=",", header = TRUE)
benchmark = read.table("SCOP_D/scores_scop_class_d.csv", sep=",", header = TRUE)
# All scores concatenated together
benchmark = read.table("tous_les_scores.csv", sep=",", header = T)

# Visualize the distribution of scores values with boxplots
boxplot(benchmark[,-c(1,2)], las=2)

# Transform the benchmarks (".", "Family", "Superfamily", "Fold") into binary numeric values:
# "Family", "Superfamily", "Fold" -> 1
#                             "." -> 0
new.Y = rep(0, dim(benchmark)[1])
new.Y[which(benchmark[, (which(colnames(benchmark) == "benchmark"))] != ".")] = 1

# Build the data.frame with benchmarks as factors (0, 1)
df.x = benchmark[, -(which(colnames(benchmark) == "benchmark"))]
df.XY = data.frame(df.x, as.factor(new.Y))

# Transform the NAs -> 0
benchmark[is.na(benchmark)] = 0

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

# Launch Generalized Linear Model
# family = binomial: because benchmark is binary 0 or 1
fit = glm(benchmark ~ ., data = train.mat, family = "binomial")
summary(fit)
# Retrieve the new scores coefficients
coefficients = coef(fit)[1,]

# Set the new global scores based on the calculated coefficients as:
# sum_scores = intersect +
#              coef1*alignment_score +
#              coef2*threading_score +
#              coef3*modeller_score +
#              coef4*secondary_structure +
#              coef5* solvent_access_score +
#              coef6*co_evolution_score
benchmark2 = benchmark
benchmark2$sum_scores = coefficients[1]+
                        coefficients[2]*benchmark2$alignment+
                        coefficients[3]*benchmark2$threading+
                        coefficients[4]*benchmark2$modeller+
                        coefficients[5]*benchmark2$secondary_structure+
                        coefficients[6]*benchmark2$solvent_access+
                        coefficients[7]*benchmark2$co_evolution
  
# Normalisation with min-max Z-Score: 0 <= score <= 1 
benchmark2$sum_scores = (benchmark2$sum_scores-min(benchmark2$sum_scores))/(max(benchmark2$sum_scores)-min(benchmark2$sum_scores))

# Renormalisation: Z-Score robuste (No giving better results)
# m <- median(benchmark2$sum_scores)
# s <- mad(benchmark2$sum_scores)
# Check doc for more infos
# The default constant = 1.4826 (approximately 1/ Φ^(-1)(3/4) = 1/qnorm(3/4)) ensures consistency,
# i.e., E[mad(X_1,…,X_n)] = σ
# for X_i distributed as N(μ, σ^2) and large n.
# Translation: The constant 1.4826 is a correction factor which makes the MAD consistent at gaussian distributions.
# benchmark2$sum_scores <- abs(benchmark2$sum_scores - m) / (s)

# Calculate the performances of the model

yobsApp=train.mat$benchmark
ppredApp=predict(fit, newdata=train.mat, type = "class")

yobsTest=test.mat$benchmark
ppredTest=predict(fit, newdata=test.mat, type = "class")

yPredApp = ppredApp
yPredTest = ppredTest

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






# ###################### SVM #####################
# 
# benchmark = read.table("tous_les_scores.csv", sep=",", header = T)
# 
# 
# benchmark = benchmark[,-c(1, 9)]
# benchmark[is.na(benchmark)] = 0
# 
# benchmark = benchmark[which(benchmark[,1] %in% c("Family", "Superfamily", "Fold")),]
# 
# # Echantillons d'apprentissage et test
# v.ind.train = sample(1:nrow(benchmark), size=2/3*nrow(benchmark), replace=FALSE)
# v.ind.test = (1:nrow(benchmark))[-v.ind.train]
# train.mat = benchmark[v.ind.train,]
# test.mat = benchmark[v.ind.test,]
# # train.mat = train.mat[,-1]
# # test.mat = test.mat[,-1]
# 
# # Construction du modèle
# result = svm(benchmark~., data=benchmark, type="C-classification", kernel="radial")
# 
# # Performances du modèle en apprentissage
# ypredApp = predict(result, train.mat[,-1])
# table(train.mat[,1], ypredApp)
# 
# # Performances du modèle en test
# ypredTest = predict(result, test.mat[,-1])
# table(test.mat[,1], ypredTest)
# 
# 
# tune.res = best.svm(benchmark~., data=benchmark, kernel="radial")
# summary(tune.res)
# plot(tune.res)
# 
# # Construction du modèle après tune
# result = svm(as.matrix(train.mat[,-10]), train.mat[,10], type="C-classification", cost=1e05, gamma=1e-04)
# 
# # Performances du modèle en apprentissage
# ypredApp = predict(result, train.mat[,-10])
# table(train.mat[,10], ypredApp)
# 
# # Performances du modèle en test
# ypredTest = predict(result, test.mat[,-10])
# table(test.mat[,10], ypredTest)