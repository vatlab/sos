## python lasso.py train.txt test.txt 5 pred.txt coef.txt
import sys
import numpy as np
from sklearn.linear_model import LassoCV
train = np.genfromtxt(sys.argv[1], delimiter = ",")
test = np.genfromtxt(sys.argv[2], delimiter = ",")
model = LassoCV(cv = int(sys.argv[3]), fit_intercept = False).fit(train[:,1:], train[:,1])
Ypred = model.predict(test[:,1:])
np.savetxt(sys.argv[4], Ypred)
np.savetxt(sys.argv[5], model.coef_)
