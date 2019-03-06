import numpy as np
import pandas as pd
import random
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.manifold import TSNE
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
random.seed(42)

 
def triple_split(filename, val_size=0.1, test_size=0.1):
    """ function to split dataset into train, validation and test sets """
    df=pd.read_csv(filename)
    train,test=train_test_split(df, test_size=test_size, shuffle=True)
    train,val=train_test_split(train, test_size=val_size, shuffle=True)
    return train,val,test


def set_prepare(df, target):
    """ function to split set into features and target """
    if ~isinstance(df,pd.DataFrame):
        df=pd.read_csv(df)
    y = df[target]
    X = df.drop([target],axis=1)
    return X,y


def ml_optimize(X,y,classifier,grid,cv=5,scoring='accuracy'):
    """ function to select hyperparameters for RandomForestClassifier """
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

 
def prot_tSNE(X,y,label='Fluor',plot=True):
    """ function to apply t-SNE dimensionality
    reduction to featurized protein set """
    tsne=TSNE(n_components=2,verbose=1,perplexity=40,random_state=42,
              n_iter=1500)
    tsne_results=tsne.fit_transform(X)
    X['x-tsne']=tsne_results[:,0]
    X['y-tsne']=tsne_results[:,1]
    df=pd.concat([X,y],axis=1)
    if plot:
        sns.set(style='white',color_codes=True)
        plt.figure()
        plt.scatter(df.loc[df[label]==0,['x-tsne']],
                    df.loc[df[label]==0,['y-tsne']],
                    marker='o',color='b',linewidth='1',alpha=0.6,
                    label='Before Transfer')
        plt.scatter(df.loc[df[label]==1,['x-tsne']],
                    df.loc[df[label]==1,['y-tsne']],
                    marker='o',color='orange',linewidth='1',alpha=0.6,
                    label='After Transfer')
        plt.xlabel('Dim 1')
        plt.ylabel('Dim 2')
        plt.title('t-SNE labeled by origin')
        plt.legend(loc='best')
        plt.show()
    else:
        return X

def rf_train_predict(X_train,y_train,X,y):
    """ takes train + test data to train a random forest
    and outputs the predictions for the test data """
    rf=RandomForestClassifier(n_estimators=100,random_state=42)
    rf.fit(X_train,y_train)
    X=X.drop(['Instability'],axis=1)
    return rf.predict(X[y==0]), rf.predict(X[y==1])
    
