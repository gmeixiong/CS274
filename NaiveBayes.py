import pandas
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import confusion_matrix
from sklearn.feature_selection import mutual_info_classif

data = pandas.read_csv("features.csv").values
for row in range(len(data)):
	for col in range(len(data[0])-1):
		if data[row][col] == 'y':
			data[row][col] = 1
		else:
			data[row][col] = 0

dataframes = data[:, [range(len(data[0])-1)]]
y = data[:, len(data[0])-1]

X = []
for i in range(len(dataframes)):
	X.append(dataframes[i][0])


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1115)
clf = BernoulliNB()


# print X_train
# print y_train
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

mutual_info = mutual_info_classif(X,y, discrete_features = True)

summ = 0
for i in range(len(mutual_info)):
	if i % 5 == 4:
		summ += mutual_info[i]
print summ

# print mutual_info.argsort()[-5:][::-1]

# print y_test
# print y_pred

# print fp
# print tn
# print fp *1.0 / (fp+tn)
