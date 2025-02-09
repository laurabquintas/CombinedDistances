
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
from matplotlib import colors as colors
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, cut_tree
import scipy.spatial as sp




def histograms(data, output_directory, method):
    bins = np.arange(0, 1.1, 0.05)

    # Get the upper triangle of the distance matrix
    upper_triangle = []
    for i in range(len(data)):
        for j in range(i+1, len(data)):
            upper_triangle.append(data[i,j])

    # Flatten the distance matrix into a 1D array
    #distances = data.flatten()

    # Plot the histogram
    plt.figure(figsize=(10,10))
    plt.hist(upper_triangle, bins=bins)
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of {method}')

    # Save the histogram
    plt.savefig(output_directory,dpi=200, transparent=True)

    # Clear the plot
    plt.clf()


def cluster_map(method, link, data, output_directory):
    data2 = sp.distance.squareform(data)
        
    #linkage
    linkage1 = linkage(data2, method=link)
    #clustermap
    cluster = sns.clustermap(data, cmap = "mako_r", cbar_pos=(1.1, .2, .03, .4), row_linkage=linkage1, col_linkage=linkage1)
    
        
    #save the clustermap
    cluster.fig.suptitle(f'Cluster Map of {method} with {link} linkage method')
    cluster.savefig(output_directory,dpi= 300, transparent=True)


def clustermap_wClinical(data, method, clinical_data, link, output_directory):
    clinical_cols = ['Overall_Tumor_Grade', 'Vital_Status', 'Receptor_Status_Patient', 'Metastatic_Dz', 'Stage']
    types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2']
    color_dict = {}


    color_dict = {}
    for (i,j) in zip(clinical_cols,types ):
        d = clinical_data[i]
        color = sns.color_palette(j, n_colors=len(d.unique()))
        lut = dict(zip(d.unique(), color))
        row_color = d.map(lut)
        color_dict[i] = row_color

    color_df = pd.DataFrame(color_dict)

    linkage1 = linkage(sp.distance.squareform(data), method=link)
    g = sns.clustermap(data, cmap='mako_r',cbar_pos=(1.1, .2, .03, .4), row_colors=color_df, row_linkage=linkage1, col_linkage=linkage1)

    g.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {link} linkage method')
    g.savefig(output_directory,dpi= 300, transparent=True)




def clustermap_wClinical_legend(data, method, clinical_data, link, output_directory):
    clinical_cols = ['PriorMalig', 'Diagnosis', 'Gender','Tx_group','VitalStatus']
    types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2']
    color_dict = {}

    for (i, j) in zip(clinical_cols, types):
        d = clinical_data[i]
        color = sns.color_palette(j, n_colors=len(d.unique()))
        lut = dict(zip(d.unique(), color))
        row_color = d.map(lut)
        color_dict[i] = row_color

    color_df = pd.DataFrame(color_dict)

    linkage1 = linkage(sp.distance.squareform(data), method=link)
    g = sns.clustermap(data, cmap='mako_r', cbar_pos=(1.1, .2, .03, .4), row_colors=color_df, 
                       row_linkage=linkage1, col_linkage=linkage1)

    # Initial position for placing the legends below the clustermap
    legend_position = (0.6, -0.05)
    
    # Add separated legend for each clinical data color below the main figure
    for (i, j) in zip(clinical_cols, types):
        d = clinical_data[i]
        color = sns.color_palette(j, n_colors=len(d.unique()))
        lut = dict(zip(d.unique(), color))
        
        # Create patches for the legend
        patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=value, markersize=10, label=key) 
                   for key, value in lut.items()]
        
        # Add legend with title
        g.fig.legend(handles=patches, loc='lower center', bbox_to_anchor=legend_position,
                     frameon=False, title=i, ncol=len(lut))
        
        # Adjust position for next legend
        legend_position = (legend_position[0], legend_position[1] - 0.05)

    g.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {link} linkage method')
    g.savefig(output_directory, dpi=300, transparent=True)

# Note: The

def clustermap_wMutations(data, method, clinical_data, link, output_directory):
    clinical_cols = ['TP53', 'FLT3', 'NPM1','NRAS','IDH2', 'NRAS','KRAS']
    types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2','Spectral', 'Set2']
    color_dict = {}

    for (i, j) in zip(clinical_cols, types):
        d = clinical_data[i]
        color = sns.color_palette(j, n_colors=len(d.unique()))
        lut = dict(zip(d.unique(), color))
        row_color = d.map(lut)
        color_dict[i] = row_color

    color_df = pd.DataFrame(color_dict)

    linkage1 = linkage(sp.distance.squareform(data), method=link)
    g = sns.clustermap(data, cmap='mako_r', cbar_pos=(1.1, .2, .03, .4), row_colors=color_df, 
                       row_linkage=linkage1, col_linkage=linkage1)

    # Initial position for placing the legends below the clustermap
    legend_position = (0.6, -0.05)
    
    # Add separated legend for each clinical data color below the main figure
    for (i, j) in zip(clinical_cols, types):
        d = clinical_data[i]
        color = sns.color_palette(j, n_colors=len(d.unique()))
        lut = dict(zip(d.unique(), color))
        
        # Create patches for the legend
        patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=value, markersize=10, label=key) 
                   for key, value in lut.items()]
        
        # Add legend with title
        g.fig.legend(handles=patches, loc='lower center', bbox_to_anchor=legend_position,
                     frameon=False, title=i, ncol=len(lut))
        
        # Adjust position for next legend
        legend_position = (legend_position[0], legend_position[1] - 0.05)

    g.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {link} linkage method')
    g.savefig(output_directory, dpi=300, transparent=True)

# Note: The

def get_cluster_labels(X, n_cluster):

    Z = linkage(sp.distance.squareform(X), method='ward')

    cluster_labels = cut_tree(Z, n_clusters=n_cluster).flatten()

    return cluster_labels
