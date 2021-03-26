


import h5py
import numpy as np
import pandas as pd

from seak.data_loaders import VariantLoader


class BurdenLoaderHDF5(VariantLoader):

    '''
    Loader class for the hdf5 files generated by export_burden.py
    '''

    def __init__(self, h5, iid, gid):

        self.iid_path = iid
        self.gid_path = gid

        with open(iid, 'r') as infile:
            iid = np.array([l.rstrip() for l in infile])

        self.iid_index = pd.Index(iid)

        with open(gid, 'r') as infile:
            gid = np.array([l.rstrip() for l in infile])

        with h5py.File(h5, 'r') as infile:
            assert (len(infile['G']) == len(gid)), "Error: length of index ({}) does not match data length ({})".format(
                len(infile['G']), len(gid))
            assert (infile['G'].shape[1] == len(
                iid)), "Error: number of individuals ({}) does not match data ({})".format(infile['G'].shape[1],
                                                                                           len(iid))

        self.gid_index = pd.Index(gid)

        self.h5_path = h5
        self.data = h5py.File(h5, 'r')['G']

        self.mask = np.arange(len(iid))

    def genotypes_by_region(self, coordinates):
        return self.data[self.gid_index.get_loc(coordinates['name']), :][np.newaxis, self.mask].T

    def genotypes_by_id(self, vids):

        if isinstance(vids, str):
            vids = np.array([vids])

        data = np.concatenate([self.data[self.gid_index.get_loc(i), :][np.newaxis, self.mask] for i in vids], axis=0).T

        return data

    def update_individuals(self, vids, exclude=False):

        if exclude:
            raise NotImplementedError

        self.mask = [self.iid_index.get_loc(i) for i in vids]

    def get_vids(self):
        return self.gid_index.values


    def get_iids(self):
        return self.iid_index.values[self.mask]
