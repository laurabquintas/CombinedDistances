import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
from matplotlib import colors as colors
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, cut_tree
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test, logrank_test
from lifelines import KaplanMeierFitter, CoxPHFitter
import scipy.spatial as sp
import io
from contextlib import redirect_stdout
import statsmodels.stats.multitest as multitest
import ast
from itertools import combinations
from sklearn.metrics import calinski_harabasz_score as chs
from sklearn.metrics import silhouette_samples, silhouette_score
import copy
import matplotlib.cm as cm
from statsmodels.nonparametric.smoothers_lowess import lowess



def silhouette_plot(method, link, data, output_directory, n_max):
    # Set the range of cluster numbers to evaluate

    Z = linkage(sp.distance.squareform(data), method=link)

    # Calculate silhouette scores for each cluster number
    silhouette_scores = []

    for n_clusters in range(2, n_max + 1):


        labels = cut_tree(Z, n_clusters=n_clusters).flatten()
        score = silhouette_score(data, labels, metric='precomputed')
        silhouette_scores.append(score)


    line = sns.scatterplot(x=range(2, n_max + 1), y=silhouette_scores)
    
    line.set(title=f'Silhouette Score of {method} with {link} linkage method')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.xticks(range(2, n_max + 1))
    plt.savefig(output_directory,dpi= 200, transparent=True)
    plt.clf()



def silhouette_plot_in_detail(method, n_clust, distance_matrix, output_directory):

    X = distance_matrix

    fig, ax1 = plt.subplots(1, 1)
    fig.set_size_inches(7, 7)

    ax1.set_xlim([-0.1, 1])

    ax1.set_ylim([0, len(X) + (n_clust + 1) * 10])

    Z = linkage(sp.distance.squareform(X), method='ward')

    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed')
    
    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='precomputed')

    y_lower = 9

    for i in range(n_clust):

        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clust)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        y_lower = y_upper + 10  

    ax1.set_title(f"The silhouette plot of {method} for {n_clust} clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    plt.savefig(output_directory, dpi=300, bbox_inches='tight')
    plt.clf()



def weighted_silhouette(n_clust, X):
    
    Z = linkage(sp.distance.squareform(X), method='ward')

    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='precomputed')

    final_silhouette = 0

    for i in range(n_clust):
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_mean = np.median(ith_cluster_silhouette_values)
        
        weight_silhouette = ith_cluster_silhouette_mean * len(ith_cluster_silhouette_values) / len(X)
        final_silhouette += weight_silhouette
    
    return final_silhouette



def get_dict_surv(methods):
    silh_surv =dict()
    for i in methods:
        silh_surv[i] = {'Clusters' : int, 'Weighted Silhouette Score': float, 'Calinski Harabasz Score': float, 'PH Assumptions' : str, 'Log Likelihood Ratio': float}
    return silh_surv



def line_plot(function, max_x, output_directory, title, x_label, y_label, color):
    line = sns.lineplot(x=range(2, max_x), y=function, marker='o', color=color)
            
    line.set(title=title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(range(2, max_x))
    plt.savefig(output_directory, dpi=300)
    plt.clf()



def show_plot(function, min_x, max_x, title, x_label, y_label, color):
    line = sns.lineplot(x=range(min_x, max_x), y=function, marker='o', color=color)
            
    line.set(title=title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(range(min_x, max_x))
    plt.show()
    plt.clf()


def prev_cluster_increase( n_clust, X, Z, patients_silhouette = None, patients_cluster = None):

    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='precomputed')

    if n_clust == 2:

        cluster_0_list = sample_silhouette_values[cluster_labels == 0]
        cluster_0 = np.median(cluster_0_list)
        cluster_0 = len(cluster_0_list) * cluster_0 / len(X)
        cluster_1_list = sample_silhouette_values[cluster_labels == 1]
        cluster_1 = np.median(cluster_1_list)
        cluster_1 = len(cluster_1_list) * cluster_1 / len(X)
        
        patients_cluster = {cluster_labels.tolist().index(0) : len(cluster_0_list), cluster_labels.tolist().index(1) : len(cluster_1_list)}
        patients_silhouette = {cluster_labels.tolist().index(0) : cluster_0, cluster_labels.tolist().index(1) : cluster_1}
        
        return patients_silhouette, patients_cluster
        
    else:
        new_patients_silhouette = dict()
        new_patients_cluster = dict()
        rate = 0
        prev_cluster = 0
        others = 0
        prev_others = 0

        for i in range(n_clust):
            cluster_list = sample_silhouette_values[cluster_labels == i]
            index = cluster_labels.tolist().index(i)
            numb = len(cluster_list)
            new_patients_cluster[index] = numb
            cluster = np.median(cluster_list)
            cluster = len(cluster_list) * cluster / len(X)
            new_patients_silhouette[index] = cluster
            if index not in patients_cluster or numb != patients_cluster[index]:
                
                rate += cluster
                
                if index in patients_silhouette:
                    prev_cluster = patients_silhouette[index]
                    for j in patients_silhouette:
                        if j != index:
                            prev_others += patients_silhouette[j]
            
            else:
                others += cluster

        new_rate = rate / (others + rate)
        
        before = prev_cluster / (prev_others + prev_cluster)
        
        if new_rate > before:
            ans = abs(new_rate - before)
        else:
            ans = -abs(before - new_rate)
    
        return ans, new_patients_silhouette, new_patients_cluster
        
            

def min_cluster_increase(n_clust :  int, X, patients_cluster : list):

    Z = linkage(sp.distance.squareform(X), method='ward')

    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='precomputed')

    if n_clust == 2:


        cluster_0_list = sample_silhouette_values[cluster_labels == 0]
        cluster_0 = np.median(cluster_0_list)
        cluster_0 = len(cluster_0_list) * cluster_0 / len(X)
        cluster_1_list = sample_silhouette_values[cluster_labels == 1]
        cluster_1 = np.median(cluster_1_list)
        cluster_1 = len(cluster_1_list) * cluster_1 / len(X)
        
        rate =min([cluster_0/cluster_1, cluster_1/cluster_0])
        
        patients_cluster = {cluster_labels.tolist().index(0) : len(cluster_0_list), cluster_labels.tolist().index(1) : len(cluster_1_list)}

        
    else:

        new_patients_cluster = dict()
        rate = list()
        others = 0
        for i in range(n_clust):
            cluster_list = sample_silhouette_values[cluster_labels == i]
            index = cluster_labels.tolist().index(i)
            numb = len(cluster_list)
            new_patients_cluster[index] = numb
            if index not in patients_cluster or numb != patients_cluster[index]:
                cluster = np.median(cluster_list)
                cluster = len(cluster_list) * cluster / len(X)
                rate.append(cluster)

            else:
                cluster = np.median(cluster_list)
                cluster = len(cluster_list) * cluster / len(X)
                others += cluster

        assert len(rate) == 2

        rate_0 = rate[0]
        rate[0] = rate[0] / (others + rate[1])
        rate[1] = rate[1] / (others + rate_0)

        #if rate[0] < 0:
            #silh_rate.append(rate[1])
        #elif rate[1] < 0:
            #silh_rate.append(rate[0])
        #else: 

        rate = min(rate)

    return rate, new_patients_cluster


def cummulative_values(list):
    cummulative = []
    cummulative.append(list[0])
    for i in range(1, len(list)):
        cummulative.append(list[i] + cummulative[i-1])
    return cummulative


def difference_dict(total_silh):

    difference = dict()

    for i in total_silh.keys():

        if i not in difference:
            difference[i] = list()

        for j in range( len(total_silh[i])-1):

            if total_silh[i][j] > total_silh[i][j+1]:
                difference[i].append(abs(total_silh[i][j]-total_silh[i][j+1]))
            else:
                difference[i].append(-abs(total_silh[i][j]-total_silh[i][j+1]))
    
    return difference

def difference(silh):
    
    difference = list()

    for i in range(len(silh)-1):

        if silh[i+1] > silh[i]:
            difference.append(abs(silh[i]-silh[i+1]))
        else:
            difference.append(-abs(silh[i]-silh[i+1]))
    
    return difference


def check_tree_clusters(tree_map, X, n_clust, n_trees):

    Z = linkage(sp.distance.squareform(X), method='ward')
    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    cluster_dict = {}
    for i in range(len(cluster_labels)):
        if cluster_labels[i] not in cluster_dict:
            cluster_dict[cluster_labels[i]] = [i+1]
        else:
            cluster_dict[cluster_labels[i]].append(i+1)

    for n in range(n_clust):
        print('Cluster ' + str(n+1) + ': \n')
        for i in range(n_trees):
            for j in tree_map[cluster_dict[n][i]].values():
                print(j)


def silhouette_except_worst(X, n_clust):

    Z =linkage(sp.distance.squareform(X), method='ward')
    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()

    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric='precomputed')

    silhouette_list = list()
    for i in range(n_clust):
        cluster_list = sample_silhouette_values[cluster_labels == i]
        cluster = np.median(cluster_list)
        cluster = len(cluster_list) * cluster / len(X)
        silhouette_list.append(cluster)
    
    silhouette_list.sort()
    silhouette_list.pop(0)

    return sum(silhouette_list)


def number_patients_new_cluster(X,n_clust):

    Z = linkage(sp.distance.squareform(X), method='ward')

    prev_cluster_labels = cut_tree(Z, n_clusters=n_clust-1).flatten()

    cluster_labels = cut_tree(Z, n_clusters=n_clust).flatten()
    
    values = list()
    prev_patients_cluster = dict()
    new_patients_cluster = dict()

    for i in range(n_clust-1):

        cluster_list = prev_cluster_labels[prev_cluster_labels == i]
        index = prev_cluster_labels.tolist().index(i)
        numb = len(cluster_list)
        prev_patients_cluster[index] = numb

    for i in range(n_clust):

        cluster_list = cluster_labels[cluster_labels == i]
        index = cluster_labels.tolist().index(i)
        numb = len(cluster_list)
        new_patients_cluster[index] = numb

        if index not in prev_patients_cluster or numb != prev_patients_cluster[index]:
            values.append(numb)

    assert len(values) == 2

    min_value = min(values)

    return min_value


def loess(X, method, x_min, x_max, output_dir=None):
    x = np.array(range(x_min, x_max+1))

    y = get_prev_inc_values(X, x_min, x_max)

    # Apply LOESS smoothing
    smoothed = lowess(y, x)

    # Plot the smoothed curve
    plt.figure()
    plt.scatter(x, y, label='Data')
    plt.plot(smoothed[:, 0], smoothed[:, 1], 'r', label='LOESS', color='darkorange')
    plt.legend()
    plt.xlabel('Number of clusters')
    plt.ylabel('Next Cluster Relative Increase (%)')
    plt.xticks(range(x_min, x_max+1))
    plt.title(f'LOESS Smoothing for {method}')
    if output_dir is None:
        plt.show()
    else:
        plt.savefig(output_dir + f'/{method}_loess.png', dpi=300, bbox_inches='tight')

def loess_simple(X, method, x_min, x_max, output_dir=None):
    x = np.array(range(x_min, x_max+1))

    y = dif_wss(X, x_max)

    # Apply LOESS smoothing
    smoothed = lowess(y, x)

    # Plot the smoothed curve
    plt.figure()
    plt.scatter(x, y, label='Data')
    plt.plot(smoothed[:, 0], smoothed[:, 1], 'r', label='LOESS', color='darkorange')
    plt.legend()
    plt.xlabel('Number of clusters')
    plt.ylabel('WSS difference')
    plt.xticks(range(x_min, x_max+1))
    plt.title(f'LOESS Smoothing for {method}')
    if output_dir is None:
        plt.show()
    else:
        plt.savefig(output_dir + f'/{method}_loess.png', dpi=300, bbox_inches='tight')

def n_clusters_loess(values, x_min, x_max, threshold=0.05):
    
    x = np.array(range(x_min, x_max+1))

    y = values

    smoothed = lowess(y, x)
    constant_indices = np.where(np.abs(smoothed[:, 1]) < threshold)[0]
    constant_x = smoothed[constant_indices[0], 0] if constant_indices.size > 0 else None

    return constant_x

def n_clusters_wss(X, min_x, max_n_clust):
    values =[]
    for n_clust in range(min_x, max_n_clust+1):
        values.append(weighted_silhouette(n_clust, X))

    index_max = values.index(max(values)) + min_x
    
    return index_max

def get_n_clusters(X, min_x = 3, max_x = 20, threshold = 0.05):

    values = get_prev_inc_values(X, max_x = max_x)

    
    loess_n = n_clusters_loess(values, min_x, max_x, threshold)
    wss_n = n_clusters_wss(X, min_x, max_x)
    out = loess_n if loess_n is not None and loess_n < wss_n else wss_n

    return out

def get_prev_inc_values(X, min_x = 3, max_x = 20):
    Z = linkage(sp.distance.squareform(X), method='ward')
    values = []
    patients_silh, patients_cluster = prev_cluster_increase(2,X,Z)
    for n_clust in range(min_x, max_x+1):
        val, patients_silh,patients_cluster = prev_cluster_increase(n_clust,X,Z,patients_silh,patients_cluster)
        values.append(val)
    
    return values
 
 
def dif_wss(X, n_clust):

    silh = list()

    for i in range(3, n_clust+2):
        silh.append(weighted_silhouette(i, X))

    dif = difference(silh)

    return dif

def get_n_clusters_simple(X, min_x = 3, max_x = 20, threshold = 0.02):

    values = dif_wss(X, max_x)
    loess_n = n_clusters_loess(values, min_x, max_x, threshold)

    wss_n = n_clusters_wss(X, min_x, max_x)

    out = loess_n if loess_n is not None and loess_n < wss_n else wss_n

    return out