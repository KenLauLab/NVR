# NVR is a python implementation of gene feature selection used in Welch et al.,2016
# Copyright (C) 2017  Bob Chen

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

import numpy as np
import scipy.spatial as sps
import networkx as nx
import time

#################################################

def parseNoise(rawArray): 
    '''
    Function returns indices that contain non-noisy genes as an integer array
    :param rawArray: numpy ndarray of data set
    :return nnGenes : numpy ndarray of non-noise gene indices
    '''
    nnGenes=[] 
    for i in range(0,(rawArray.shape[1])): #Checks all genes
        count0=np.asarray(np.where(rawArray[:,i]==0))
        count1=np.asarray(np.where(rawArray[:,i]==1))
        if ((count1.shape[1]+count0.shape[1])<rawArray.shape[0]):
            nnGenes=np.append(nnGenes,i)
    return nnGenes.astype(int)

#################################################

def mkIndexedArr(unparsedArr,nnGeneIndex):
    '''
    Function returns ndarray of gene expression values from ndarray of unmodified data set and ndarray of gene indices
    :param unparsedArr: numpy ndarray of unmodified data set
    :param nnGeneIndex: numpy ndarray of indices of genes of interest
    :return nnGeneExp:  numpy ndarray of data set including all cells and expression data only of genes of interest 
    '''
    nnGeneExp=np.zeros([unparsedArr.shape[0],nnGeneIndex.shape[0]],dtype=float) #initializes new array for output
    for i in range(0,nnGeneIndex.shape[0]): #fill in the new array based on indices of interest
        nnGeneExp[:,i]=unparsedArr[:,nnGeneIndex[i]]
    return nnGeneExp

#################################################

def pwArcsinh(inputArray,constant): 
    '''
    Function returns an ndarray by performing a pointwise inverse hyperbolic sine transformation on the input ndarray 
    :param inputArray: ndarray of gene expression data in raw count or RPKM format 
    :param contstant: some constant to normalize for total read count number
    :return transformed: ndarray of transformed inputArray 
    '''
    #this function assumes rows are cells and columns are genes in the input array
    countsum=np.sum(inputArray,axis=1) #calculate the total number of counts per cell
    holder=np.zeros_like(inputArray,dtype=float)
    transformed=np.zeros_like(inputArray,dtype=float) #initialize the output array
    print ("Completion:")
    for i in range(0,inputArray.shape[0]):
        #print (((i/inputArray.shape[0])*100),end='\r')#progress meter 
        holder[i,:]=(inputArray[i,:]/countsum[i])*constant #divide each genes counts by total number of counts in cell
        transformed[i,:]=np.arcsinh((holder[i,:]/constant*1000)) #do arcsinh transform for each element of matrix
    return transformed

#################################################

def adaptive_knn_graph(traj_dist,k): 
    '''
    Function returns an adjacency matrix based on the euclidean distance of each gene's expression across all cells
    :param traj_dist: ndarray of euclidean distances
    :param k: int of the minimum number of connections needed for a minimally connected graph, calculated by min_conn_k()
    :return adj_mat: ndarray representing the calculated adjacency matrix
    '''
    adj_mat = np.zeros_like(traj_dist,dtype=float) #initialize output matrix
    knn=(np.transpose(np.argsort(traj_dist,0))+1) #eturns the indices that would sort an array, transposed for consistent formatting
    for i in range(0,(traj_dist.shape[0])): #go through the whole distance matrix and set the corresponding adjacencies
        adj_mat[i,knn[i,range(1,k+1)]-1]=traj_dist[i,knn[i,range(1,k+1)]-1]
    return adj_mat

#################################################

def min_conn_k(traj_exp):
    '''
    Function returns the minimum number of connections, k, that are required to form a fully connected graph based gene expression data
    :param traj_exp: ndarray representing gene expression
    :return k: int of the minimum number of connections needed for a minimally connected graph
    '''
    traj_dist=sps.distance.squareform(sps.distance.pdist(traj_exp)) #this is the gene expression euclidean distance calculation
    conn_comp=2 #just a starting value for the number of connected components in the graph
    k=0 #just a starting value for the number of neighbors required in the graph
    while (conn_comp>1):
        k=k+1 #each time this loops, increase the number of neighbors by 1 until we get 1 graph component(indicating a fully connected graph)
        adj_mat=adaptive_knn_graph(traj_dist,k) #uses adaptive_knn_graph to generate an adjacency matrix
        traj_graph = nx.Graph(adj_mat) #uses that adjacency matrix to make a graph
        conn_comp = nx.number_connected_components(traj_graph) #count the number of connected components in that graph
    return k

#################################################

def dev_ij(i,j,traj_exp,adj_mat): #takes adjacency matrix and raw expression data
    '''
    Function returns the sum of squared differences in the expression of a gene, per cell
    :param i: int representing the index of the cell being processed 
    :param j: int representing the index of the gene being processed
    :param traj_exp: ndarray representing gene expression
    :param adj_mat: ndarray representing the calculated adjacency matrix
    :return: float representing the sum of squared differences
    '''
    t=np.asmatrix(traj_exp) 
    wh=np.where(adj_mat[i]>0) #looks for indices where the adjacency is greater than 0
    return np.sum(np.square(t[i,j]-t[wh,j],dtype=float),dtype=float) #performs sum of squared difference calculation

#################################################

def selection_val(traj_exp,adj_mat):
    '''
    Function returns an ndarray of ratios calculated by dividing the summed neighborhood variances by the global variance
    :param traj_exp: ndarray representing gene expression
    :param adj_mat: ndarray representing the calculated adjacency matrix
    :return val: ndarray representing the ratio of summed neighborhood variances by the global variance 
    '''
    r = traj_exp.shape[0] #keep track of the rows
    c = traj_exp.shape[1] #keep track of the columns
    k = np.sum(adj_mat[0]>0,dtype=float) #k calculation based off of adjacency matrix
    dev=np.zeros_like(traj_exp,dtype=float) #initialize matrix to store dev_ij values
    val=np.zeros(traj_exp.shape[1],dtype=float) #initialize val matrix to store variance values
    print ("Start global variance calculation")
    ivar=np.var(traj_exp,axis=0,ddof=1) #calculate variance in traj_exp, this is the global variance matrix
    print ("Finished global variance calculation")
    print ("Start neighborhood variance calculation") 
    print ("Completion:")
    for i in range(0,r): #start of dev_ij calculation loop
        #print (((i/r)*100),end='\r') #progress meter
        for j in range(0,c):
            dev[i,j] = dev_ij(i,j,traj_exp,adj_mat) #this is the part that takes the longest
    print ("Finished neighborhood variance calculation")
    rowsumdev=np.sum(dev,axis=0) #axis-wise sum of deviations calculated between i and j
    print ("Start global to neighborhood variance ratio calculation")
    for i in range(0,c): #set values to variance/deviation calculation in loop
        if rowsumdev[i]!=0: #pretty much just throw out anything that has devsum=0, when considering deviation sums of 0, we could get a divide by zero error
            val[i] = ((ivar[i])/(rowsumdev[i]/(r*k-1))) #variance calculation done here and written into val matrix
    print ("Finished global to neighborhood variance ratio calculation")
    return val

#################################################

def select_genes(embedding):
    '''
    Wrapper function to perform neighborhood variance based feature selection
    :param embedding: ndarray representing gene expression, consists of cells as rows and genes as columns
    :return genes_mat: ndarray of selected genes
    '''
    print ("Start min_conn_k")
    start=time.time() #just start a timer to see how long this whole thing takes in seconds
    k = min_conn_k(embedding) #run the min_conn_k to find k
    print (k, "connections needed") #prints the k for reference
    print ("Finished min_conn_k ")
    r = embedding.shape[0]
    c = embedding.shape[1]
    print ("Start traj_dist")
    traj_dist = sps.distance.squareform(sps.distance.pdist(embedding)) #distance calculation to use in making the adjacency matrix
    print ("Finished traj_dist")
    print ("Start adaptive_knn_graph")
    adj_mat = adaptive_knn_graph(traj_dist,k) #use traj_dist to make adjacency matrix
    print ("Finished adaptive_knn_graph")
    sel_vals = selection_val(embedding,adj_mat) #calculate selection values, or the ratios of global variance to neighborhood variance
    print ("Finished selection_val") 
    genes = np.where(sel_vals > 1) #decision made here, if neighborhood variance is lower than global, select that gene
    print ("Finished gene selection in",time.time()-start, "seconds") 
    #genes_mat = np.asmatrix(genes) #output of genes as a matrix
    print ("done")
    return np.squeeze( genes)
