#!/usr/bin/env python

"""asmd.py: ACSPO SPT Mask Differences
   -----------------------------------

   This program helps visualize questionable areas within the SPT mask
   Questionable areas are defined as pixels where the ACSPO mask is labeled
   non-ocean while the SPT mask has labeled the pixels ocean.  These areas
   are visualized by running a connected component algorithm on the questionable
   areas and created images centered around the areas.

   overlay description:
   0: spt_mask == clear & acspo_mask == clear
   1: spt_mask == clear & acspo_mask == cloud
   2: spt_mask == clear & acspo_mask == front
   3: spt_mask == cloud & acspo_mask == clear
   4: spt_mask == cloud & acspo_mask == cloud
   5: spt_mask == cloud & acspo_mask == front"""


import os
import sys

import netCDF4
from skimage import measure
from scipy import ndimage
import numpy as np
from matplotlib import pyplot as plt
"""
Main steps:

1) open nc file
	a) get acspo
	b) get spt
	c) get albedo
	d) sst
2) combine acspo and spt
3) check albedo for values
4) compute gradient
5) cluster combined acspo and spt
6) set colors for ocean, front, restored, and questionable
7) display albedo,(if available), sst, grdient, overlay for each cluster with cluster in center
"""

class Granule(object):
    """A granule contains the following layers:

    Attributes:
        sst: Sea surface temperature.
        ACSPO: ACSPO cloud mask.
        albedo: the albedo band
        SPT_mask: the updated ACSPO clou dmask with fronts
        gradient: gradient magnitude of the sst
        overlay: combination of ACSPO and SPT_mask
    """

    def __init__(self, path):
        """Return a granule object with all layers"""
        if not os.path.exists(path):
        	raise RuntimeError, "File %s does not exist" % path
        self.path = path
        

    def _read_var(self, variable):
    	cdf = netCDF4.Dataset(self.path)
    	data = np.squeeze(cdf[variable][:])
    	return data

    def	_get_acspo_cloud_mask(self, acspo_flags):
    	mask_bits = 192
    	mask = np.bitwise_and(acspo_flags,mask_bits).astype(bool)
    	return mask

    def _compute_overlay(self, acspo_mask, spt_mask):
        overlay = np.zeros(acspo_mask.shape).astype(np.uint8)
        overlay[(acspo_mask==0) & (spt_mask==0)] = 0 
        overlay[(acspo_mask==0) & (spt_mask==3)] = 1
        overlay[(acspo_mask==0) & (spt_mask==4)] = 2
        overlay[(acspo_mask==1) & (spt_mask==0)] = 3 
        overlay[(acspo_mask==1) & (spt_mask==3)] = 4
        overlay[(acspo_mask==1) & (spt_mask==4)] = 5 
    	return overlay

    def _compute_gradient(self, sst):
    	dX,dY = np.gradient(sst)
    	return np.sqrt(dX**2 + dY**2)

    def get_layers(self):
    	"""open the cdf file and get the necessary variable
    	   return a dictionary with the layer names as keys and 
    	   matrices as value"""

    	layers = {}
    	layers['albedo'] = self._read_var('albedo_chM10')
    	layers['sst'] = self._read_var('sst_regression')
    	spt_mask = self._read_var('spt_mask')
    	layers['albedo'] = self._read_var('albedo_chM10')
    	
    	acspo_flags = self._read_var('acspo_mask')
    	acspo_mask =  self._get_acspo_cloud_mask(acspo_flags)

    	layers['overlay'] = self._compute_overlay(acspo_mask,spt_mask)
    	layers['gradient'] = self._compute_gradient(layers['sst'])
        layers['questionable'] = ((layers['overlay']==3) | (layers['overlay']== 5))
    	return layers





if len(sys.argv) < 2:
	print "usage: python asmd.py <granule_name>"
	sys.exit()

path = sys.argv[1]
granule = Granule(path)
layers = granule.get_layers()

#compute connected components
labels, n_lbls = measure.label(layers['questionable'],background=0, return_num=True)

#compute centers of clusters and left right up and down of cluster
centers = {}
inds_left = {}
inds_right = {}
inds_down = {}
inds_up = {}
for lbl in range(n_lbls):
    mask = labels == lbl
    centers[lbl] = ndimage.measurements.center_of_mass(labels,mask)
    inds_row, inds_col =  np.where(mask)
    inds_left[lbl] = np.min(inds_col)
    inds_right[lbl] =np.max(inds_col)
    inds_up[lbl] = np.min(inds_row)
    inds_down[lbl] = np.max(inds_row)