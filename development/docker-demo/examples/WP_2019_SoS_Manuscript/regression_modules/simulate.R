## Rscript simulate.R seed=1 N="c(40,200)" b="c(3, 1.5, 0, 0, 2, 0, 0, 0)" rstd=3 oftrain="'train.txt'" oftest="'test.txt'"
eval(parse(text=commandArgs(T)))
set.seed(seed)
p = length(b)
X = MASS::mvrnorm(n = sum(N), rep(0, p), 0.5^abs(outer(1:p, 1:p, FUN = "-")))
Y = X %*% b + rnorm(sum(N), mean = 0, sd = rstd)
Xtrain = X[1:N[1],]; Xtest = X[(N[1]+1):sum(N),]
Ytrain = Y[1:N[1]]; Ytest = Y[(N[1]+1):sum(N)]
write.table(cbind(Ytrain, Xtrain), oftrain, row.names = F, col.names = F, sep = ',')
write.table(cbind(Ytest, Xtest), oftest, row.names = F, col.names = F, sep = ',')
