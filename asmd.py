#!/usr/bin/env python

"""asmd.py: ACSPO SPT Mask Differences
   -----------------------------------

   This program helps visualize questionable areas within the SPT mask
   Questionable areas are defined as pixels where the ACSPO mask is labeled
   non-ocean while the SPT mask has labeled the pixels ocean.  These areas
   are visualized by running a connected component algorithm on the questionable
   areas and created images centered around the areas."""


import os
import sys

import netCDF4
from skimage import measure
import numpy as np
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
    	mask = np.bitwise_and(acspo_flags,mask_bits)
    	return mask

    def _compute_overlay(self, acspo_mask, spt_mask):
    	return 2*spt_mask + acspo_mask

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

    	return layers





if len(sys.argv) < 2:
	print "usage: python asmd.py <granule_name>"
	sys.exit()

path = sys.argv[1]
granule = Granule(path)
layers = granule.get_layers()



