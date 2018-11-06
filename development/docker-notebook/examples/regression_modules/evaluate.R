## Rscript evalute.R b="c(3, 1.5, 0, 0, 2, 0, 0, 0)" test="'test.txt'" fpred="'pred.txt'" fcoef="'coef.txt'" output="'score.txt'"
eval(parse(text=commandArgs(T)))
Ytruth = as.matrix(read.csv(test, header = F)[,-1]) %*% b
Ypred = scan(fpred)
prediction_mse = mean((Ytruth - Ypred)^2)
betahat = scan(fcoef)
estimation_mse = mean((betahat - b) ^ 2)
cat(paste(prediction_mse, estimation_mse), file = output)
