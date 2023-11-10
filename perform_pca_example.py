# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 08:05:11 2023

@author: mast527
"""

import pandas as pd
import numpy as np
from sklearn.datasets import make_blobs
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer




#pd.options.plotting.backend = "plotly"
pd.options.plotting.backend = "matplotlib"


import matplotlib.pyplot as plt
def simulated_data(num_samples=1000, num_features=50, centers=3):
    # Set a random seed for reproducibility
    np.random.seed(0)
    
    # Generate synthetic blobs of data using make_blobs
    num_samples = num_samples  # Number of data points
    num_features = num_features   # Number of features
    centers = centers        # Number of clusters
    
    data, labels = make_blobs(n_samples=num_samples, n_features=num_features, centers=centers, cluster_std=10)
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=[f"feature_{i+1}" for i in range(num_features)])
    
    # Add cluster labels to the DataFrame
    df['Cluster'] = labels
    
    peptides = [f'sample_{i}' for i in range(len(df))]
    df['sample'] = peptides
    df = df.set_index('sample')
    
    return df

def impute_data(X, method='median'):
    
    if method == 'mean':
        imputer = SimpleImputer(strategy='mean')
    
    if method == 'median':
        imputer = SimpleImputer(strategy='median')

    try:
        data_imputed = imputer.fit_transform(X)
    except Exception as err:
        print(err)
        print('_impute_data_failed()')
        raise Exception(err)
    return data_imputed

def perform_pca(df, feature_columns=None, category_columns:list=None, impute=False):
    '''
    
    '''
    df = df.copy()
    
    

    # Standardize the data (recommended for PCA)
    scaler = StandardScaler()
    if feature_columns is not None:
        if impute:
            X = impute_data(df[feature_columns].to_numpy())
            scaled_data = scaler.fit_transform(X)
        else:
            scaled_data = scaler.fit_transform(df[feature_columns])
    else:
        if impute:
            X = impute_data(df)
            scaled_data = scaler.fit_transform(X)
        else:
            scaled_data = scaler.fit_transform(df)
    
    # Perform PCA with the desired number of components (e.g., 2)
    num_components = 3
    pca = PCA(n_components=num_components)
    principal_components = pca.fit_transform(scaled_data)
    
    # Create a new DataFrame with the principal components
    pca_df = pd.DataFrame(data=principal_components, columns=[f"PC_{i+1}" for i in range(num_components)], index=df.index.to_list())
    
    
        
    # Get the loadings (eigenvectors) from the PCA model
    loadings = pd.DataFrame(pca.components_.T, columns=pca_df.columns, index=feature_columns)
    if category_columns is not None:
        for col in category_columns:
            pca_df[f'{col}'] = df[col].to_list() 
    return pca_df, loadings




def simulated_data2(num_samples=50, num_features=5000,centers=3):
    # Set a random seed for reproducibility
    np.random.seed(0)
    
    # Generate synthetic blobs of data using make_blobs
    num_samples = num_samples  # Number of data points
    num_features = num_features   # Number of features
    centers = centers        # Number of clusters
    
    data, labels = make_blobs(n_samples=num_samples, n_features=num_features, centers=centers, cluster_std=10)
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=[f"feature_{i+1}" for i in range(num_features)])
    
    # Add cluster labels to the DataFrame
    df['Cluster'] = labels
    
    peptides = [f'sample_{i}' for i in range(len(df))]
    df['samples'] = peptides
    df = df.set_index('samples')
    
    return df



def perform_pca2(df, feature_columns=None, category_columns:list=None):
    
    '''
    performs PCA on the samples instead of the features
    assume that features represent peptides, and samples represent experimental samples like a subject
        
    '''
    df = df.copy()
    
    # Standardize the data (recommended for PCA)
    scaler = StandardScaler()
    if feature_columns is not None:
        scaled_data = scaler.fit_transform(df[feature_columns].T)
    else:
        scaled_data = scaler.fit_transform(df)
    
    # Perform PCA with the desired number of components (e.g., 2)
    num_components = 3
    pca = PCA(n_components=num_components)
    principal_components = pca.fit_transform(scaled_data)
    
    if feature_columns is not None:
        # Create a new DataFrame with the principal components
        pca_df = pd.DataFrame(data=principal_components, columns=[f"PC_{i+1}" for i in range(num_components)], index=feature_columns)
    else:
        pca_df = pd.DataFrame(data=principal_components, columns=[f"PC_{i+1}" for i in range(num_components)])
    
        
    # Get the loadings (eigenvectors) from the PCA model
    loadings = pd.DataFrame(pca.components_.T, columns=pca_df.columns, index=df.index.to_list())
    if category_columns is not None:
        for col in category_columns:
            pca_df[f'{col}'] = df[col].to_list() 
    return pca_df, loadings

if __name__ == '__main__':
    #def PCA_method1():
    # #this simulates PCA on peptide columns
    num_samples = 5000 # Number of data points
    num_features = 500   # Number of features
    centers = 3        # Number of clusters

    d1 = simulated_data(num_samples=num_samples, num_features=num_features,centers=centers)
    feature_cols1 = d1.columns.to_list()[:num_features]
    x_pca1, x_loading1 = perform_pca(d1, feature_columns=feature_cols1, category_columns=['Cluster'])
    x_pca2, x_loading2 = perform_pca2(d1, feature_columns=feature_cols1) 

    plt.figure()
    x_pca1.plot.scatter(x='PC_1', y='PC_2', c='Cluster',cmap='viridis', title='PCA simulate_data() Method 1, perfrom_pca() Method 1')
    x_loading1.plot.scatter(x='PC_1', y='PC_2', title='Loading simulate_data() Method 1, perfrom_pca() Method 1')


    plt.figure()
    x_pca2.plot.scatter(x='PC_1', y='PC_2', title='PCA simulate_data() Method 1, perfrom_pca() Method 2')
    x_loading2.plot.scatter(x='PC_1', y='PC_2', title='Loading simulate_data() Method 1, perfrom_pca() Method 2')



    #################################################################################

    #def PCA_method2():
    num_samples = 50  # Number of data points
    num_features = 5000   # Number of features
    centers = 3        # Number of clusters

    d2 = simulated_data2(num_samples=num_samples, num_features=num_features,centers=centers)
    feature_cols2 = d2.columns.to_list()[:num_features]
    x_pca1, x_loading1 = perform_pca(d2, feature_columns=feature_cols2, category_columns=['Cluster']) 
    x_pca2, x_loading2 = perform_pca2(d2, feature_columns=feature_cols2) 
    plt.figure()

    plt.figure()
    x_pca1.plot.scatter(x='PC_1', y='PC_2', c='Cluster',cmap='viridis', title='PCA simulate_data() Method 2, perfrom_pca() Method 1')
    x_loading1.plot.scatter(x='PC_1', y='PC_2', title='Loading simulate_data() Method 2, perfrom_pca() Method 1')


    plt.figure()
    x_pca2.plot.scatter(x='PC_1', y='PC_2', title='PCA simulate_data() Method 2, perfrom_pca() Method 2')
    x_loading2.plot.scatter(x='PC_1', y='PC_2', title='Loading simulate_data() Method 2, perfrom_pca() Method 2')



