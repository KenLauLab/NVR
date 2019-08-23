# NVR is a python implementation of gene feature selection used in Welch et al.,2016
# Copyright (C) 2019  Bob Chen

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.from __future__ import division,print_function

from __future__ import division,print_function
import numpy as np
import scanpy as sc
import pandas as pd
import igraph
import time
import warnings

#################################################

def arcsinh_transform(AnnData, cofactor = 1000):
    AnnData.X = np.arcsinh(AnnData.X*cofactor,dtype='float')

#################################################

def min_conn_knn(AnnData):
    conn_comp = 0
    k = 0
    while (conn_comp!=1):
        k = k+1
        sc.pp.neighbors(AnnData,k,use_rep = 'X')
        adata_graph = igraph.Graph.Adjacency(AnnData.uns['neighbors']['distances'].toarray().astype(bool).astype(int).tolist(),mode=1) #sc.Neighbors(AnnData).to_igraph()
        conn_comp = len(list(adata_graph.components()))
        print(k,'-neighbor(s) results in',conn_comp,'component(s)')
    print('Using minimally connected KNN graph with',k,'neighbors')
    adata_graph = igraph.Graph.Adjacency(AnnData.uns['neighbors']['distances'].toarray().astype(bool).astype(int).tolist(),mode=0) #sc.Neighbors(AnnData).to_igraph()
    return(AnnData, adata_graph)
        
#################################################

def neighborhood_MSE(AnnData, graph_in):
    AnnData.obsm['Neighborhood_MSE'] = np.zeros_like(AnnData.X) #initialize neighborhood variance observation
    for i in range(graph_in.vcount()):
        neighborhood_i_values = AnnData.X[graph_in.neighborhood(i,mode='out')]
        AnnData.obsm['Neighborhood_MSE'][i] = np.sum(np.square((neighborhood_i_values[0]-neighborhood_i_values[1:])),axis=0) #sum sq diff per cell neighborhood, halve it since its only the upper triangle
    return AnnData

#################################################

def nvr_feature_select(AnnData, suppress_warnings = True):
    if suppress_warnings:
        warnings.filterwarnings("ignore")
    start=time.time() #just start a timer to see how long this whole thing takes in seconds
    AnnData,agraph = min_conn_knn(AnnData)
    print("Calculating neighborhood MSE")
    AnnData = neighborhood_MSE(AnnData,agraph)
    print("Calculating global/neighborhood variance ratio")
    AnnData.uns['Global_Variance'] = np.var(AnnData.X,axis=0,ddof=1)
    AnnData.uns['Summed_MSE'] = np.sum(AnnData.obsm['Neighborhood_MSE'],axis=0)
    AnnData.uns['Edges_per_neighborhood'] = AnnData.uns['neighbors']['distances'].toarray().astype(bool)[0].sum()
    AnnData.uns['NVR_genes'] = (np.nan_to_num(np.divide(AnnData.uns['Global_Variance'],(AnnData.uns['Summed_MSE']/(agraph.vcount()*AnnData.uns['Edges_per_neighborhood']-1))))>1)
    print("Feature selection complete,",int(time.time()-start), "seconds elapsed")
    return(AnnData)

#################################################

def subsample(partitions,dataset,seed):
    '''
    Function to generate randomly sampled datasets with replacement. This is in the context of cells in the native dataset
    which are the rows of the matrix
    :param partitions: int designating the number of evenly spaced sample sizes to randomly select from the native dataset
    :param dataset: DataFrame of the native dataset compatible with the suffle function
    :param seed: pseudorandom seed, compatible with the replicate wrapper since it adds the index to the seed
    :return subOut: dictionary of the randomly sampled datasets, keys are the number of cells
    '''
    parts=np.arange(dataset.shape[0]/partitions,dataset.shape[0],dataset.shape[0]/partitions).astype(int)    
    subOut={}
    for i in range(parts.shape[0]):
        subOut["{0}cells".format(parts[i])]=np.asarray(shuffle(dataset,random_state=seed))[0:parts[i],:]
    return subOut

#################################################

def subsampleReplicates(repNumber,partitions,dataset,seed):
    '''
    Wrapper function that generates replicate datasets using the subsampling function.
    :param repNumber: int number of replicates to generate based on the parameters given
    :param partitions: int designating the number of evenly spaced sample sizes to randomly select from the native dataset
    :param dataset: DataFrame of the native dataset compatible with the suffle function
    :param seed: pseudorandom seed, compatible with the replicate wrapper since it adds the index to the seed
    :return repOut: nested dictionary of the randomly sampled datasets, keys are the replicate number
    '''
    repOut={}
    for i in range(repNumber):
        repOut["replicate{0}".format(i)]=subsample(partitions,dataset,seed+i)
    return repOut

#################################################

def dictToFile(dictionary,replicateKey,outFileName):
    '''
    Function to write dictionary data, from subsampleReplicates, to file an hdf5 file. 
    :param dictionary: nested dictionary returned by subsampleReplicates
    :param replicateKey: string designating the replicate written to file
    :param outFileName: string defining the hdf5 filename
    '''
    replicateToFile=h5py.File(outFileName,"w")
    for i in range(len(dictionary[replicateKey])):
        replicateToFile.create_dataset("{}".format(dictionary[replicateKey].keys()[i])\
                                    ,data=dictionary[replicateKey].values()[i]\
                                    ,compression="gzip")
    replicateToFile.close()