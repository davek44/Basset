<img src="docs/basset_image.png" width="200">

# Basset
#### Deep convolutional neural networks for DNA sequence analysis.

Basset provides researchers with tools to:

1. Train deep convolutional neural networks to learn highly accurate models of DNA sequence activity such as accessibility (via DNaseI-seq or ATAC-seq), protein binding (via ChIP-seq), and chromatin state.
2. Interpret the principles learned by the model.

---------------------------------------------------------------------------------------------------
### Installation

Basset has a few dependencies because it uses both Torch7 and Python and takes advantage of a variety of packages available for both.

First, I recommend installing Torch7 from [here](http://torch.ch/docs/getting-started.html). If you plan on training models on a GPU, make sure that you have CUDA installed and Torch should find it.

For the Python dependencies, I highly recommend the [Anaconda distribution](https://www.continuum.io/downloads). The only library missing is pysam, which you can install through Anaconda or manually from [here](https://code.google.com/p/pysam/). If you don't want to use Anaconda, check out the full list of dependencies [here](docs/requirements.md).

Basset relies on the environmental variable BASSETDIR to orient itself. In your startup script (e.g. .bashrc), write
```
    export BASSETDIR=the/dir/where/basset/is/installed
```

To make the code available for use in any directory, also write
```
    export PATH=$BASSETDIR/src:$PATH
    export LUA_PATH=$BASSETDIR/src/?.lua;;LUA_PATH
    export PYTHONPATH=$BASSETDIR/src:$PYTHONPATH
```

To download and install the remaining dependencies, run
```
    ./install_dependencies.py
```

Alternatively, Dr. Lee Zamparo generously [volunteered his Docker image](https://hub.docker.com/r/lzamparo/basset/).

To download and install additional useful data, like my best pre-trained model and public datasets, run
```
    ./install_data.py
```

---------------------------------------------------------------------------------------------------
### Documentation

Basset is under active development, so don't hesitate to ask for clarifications or additional features, documentation, or tutorials.

- [File specifications](docs/file_specs.md)
  - [BED](docs/file_specs.md#bed)
  - [Table](docs/file_specs.md#table)
  - [HDF5](docs/file_specs.md#hdf5)
  - [Model](docs/file_specs.md#model)
- [Preprocess](docs/preprocess.md)
  - [preprocess_features.py](docs/preprocess.md#preprocess_features.py)
  - [seq_hdf5.py](docs/preprocess.md#seq_hdf.py)
  - [basset_sample.py](docs/preprocess.md#basset_sample.py)
- [Learning](docs/learning.md)
  - [basset_train.lua](docs/learning.md#train)
  - [basset_test.lua](docs/learning.md#test)
  - [basset_predict.lua](docs/learning.md#predict)
- [Visualization](docs/visualization.md)
  - [basset_motifs.py](docs/visualization.md#motifs)
  - [basset_motifs_infl.py](docs/visualization.md#infl)
  - [basset_sat.py](docs/visualization.md#sat)
  - [basset_sat_vcf.py](docs/visualization.md#sat_vcf)
  - [basset_sad.py](docs/visualization.md#sad)

---------------------------------------------------------------------------------------------------
### Tutorials

These are a work in progress, so forgive incompleteness for the moment. If there's a task that you're interested in that I haven't included, feel free to post it as an Issue at the top.

- Preprocess
  - [Prepare the ENCODE and Epigenomics Roadmap compendium from scratch.](tutorials/prepare_compendium.ipynb)
  - [Prepare new dataset(s) by adding to a compendium.](tutorials/new_data_many.ipynb)
  - [Prepare new dataset(s) in isolation.](tutorials/new_data_iso.ipynb)
- Train
  - [Train a model.](tutorials/train.md)
- Test
  - [Test a trained model.](tutorials/test.ipynb)
- Visualization
  - [Study the motifs learned by the model.](tutorials/motifs.ipynb)
  - [Execute an in silico saturated mutagenesis](tutorials/sat_mut.ipynb)
  - [Compute SNP Accessibility Difference profiles.](tutorials/sad.ipynb)
