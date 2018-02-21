from setuptools import setup

setup( name='NVR',
       version='0.0.2',
       description='python implementation of gene feature selection used in Welch et al.,2016',
       url='https://github.com/KenLauLab/NVR',
       author='Bob Chen',
       author_email='bob.chen@vanderbilt.edu',
       packages=['nvr'],
       install_requires=[
           'numpy>=1.11.0',
           'pandas>=0.17.1',
           'networkx',
           'scipy>=0.17.0'
       ],
     )

