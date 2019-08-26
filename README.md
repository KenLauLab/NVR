# **NVR (neighborhood variance ratio)**

Python implementation of NVR (neighborhood variance ratio) gene selection to select genes with local and monotonic variation [(Welch et al., 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27215581). The selected genes possess specific expression patterns over the entire data space amenable to trajectory analysis. Support for python 2 is now deprecated with version 0.1.0 scanpy integration, as scanpy requires python 3.

### Installation for Mac or Linux

0. 

NVR requires igraph and scanpy for operation, so please see our [instructions](https://github.com/KenLauLab/pCreode) for the installation of these packages.

There are two ways to install NVR with Mac/Linux operating systems.

1.
```python
git clone git://github.com/KenLauLab/NVR
cd NVR
pip install .
```

2.
```python
pip install nvr
```

### Installation for Windows (in Anaconda3 environment)
In your anaconda terminal, navigate to a directory where you want to install NVR (e.g. d:\src)
```
cd d:\src
git clone git://github.com/KenLauLab/NVR
cd NVR
pip install .
```

### Tutorial

Please see the NVR running .ipynb [tutorial](https://github.com/KenLauLab/NVR/blob/master/notebooks/NVR_tutorial.ipynb) as well as its associated [dataset](https://github.com/bobchen1701/NVR/tree/master/data/s1_counts.h5ad)
