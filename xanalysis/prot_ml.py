import numpy as np
import pandas as pd
import random
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
random.seed(42)

# function to split dataset into train, validation and test sets
def triple_split(filename, val_size=0.1, test_size=0.1):
    df=pd.read_csv(filename)
    train,test=train_test_split(df, test_size=test_size, shuffle=True)
    train,val=train_test_split(train, test_size=val_size, shuffle=True)
    return train,val,test

# function to split set into features and target
def set_prepare(df, target):
    y = df[target]
    X = df.drop([target],axis=1)
    return X,y

# function to select hyperparameters for RandomForestClassifier
def ml_optimize(X,y,classifier,grid,cv=5,scoring='accuracy'):
    gsearch=GridSearchCV(estimator=classifier(random_state=42),
                         param_grid=grid,scoring=scoring,
                         iid=False,cv=cv)
    gsearch.fit(X,y)
    print("Grid scores on development set:")
    means = gsearch.cv_results_['mean_test_score']
    stds = gsearch.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, gsearch.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print(gsearch.best_params_, gsearch.best_score_)
    return gsearch.best_estimator_
