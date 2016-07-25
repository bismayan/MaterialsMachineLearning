from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

def clust(dt,n_clust=10,ran_state=42):
    kmeans=KMeans(n_clusters=n_clust,random_state=ran_state)
    dt_arr=dt.values
    labels=kmeans.fit_predict(dt_arr)
    return labels



