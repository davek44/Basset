### Basset
###### Deep convolutional neural networks for DNA sequence analysis.
--------------------------------------------------------------------------------
### File specifications

<a name="bed"/>
###### BED

- chrom
- start
- end
- name (unused)
- score (unused)
- strand
- accessibilities : comma-separated list of integer indexes

--------------------------------------------------------------------------------
<a name="table"/>
###### Table

- tab-separated
- row indexes are `chrom:start:end:strand` corresponding to [BED](#bed).
- columns are sample names.
- entries are 0/1 accessibility/binding/activity.

--------------------------------------------------------------------------------
<a name="hdf5"/>
###### HDF5

- https://www.hdfgroup.org/HDF5/
- Basset will look for datasets named:
  - train_in
  - train_out
  - valid_in
  - valid_out
  - test_in
  - test_out
  - test_headers

--------------------------------------------------------------------------------
<a name="model"/>
###### Model

- [serialized Torch7 model](https://github.com/torch/torch7/blob/master/doc/serialization.md)
- [basset_train.lua](learning.md#train) saves a check point after every epoch and saves the best validation loss.
