## Rscript ridge.R train="'train.txt'" test="'test.txt'" nfolds=5 ofpred="'pred.txt'" ofcoef="'coef.txt'"
eval(parse(text=commandArgs(T)))
train = read.csv(train, header = F)
test = read.csv(test, header = F)
model = glmnet::cv.glmnet(as.matrix(train[,-1]), train[,1], family = "gaussian", alpha = 0, nfolds = nfolds, intercept = F)
betahat = as.vector(coef(model, s = "lambda.min")[-1])
Ypred = predict(model, as.matrix(test[,-1]), s = "lambda.min")
write.table(Ypred, ofpred, row.names = F, col.names = F, sep = ',')
write.table(betahat, ofcoef, row.names = F, col.names = F, sep = ',')
