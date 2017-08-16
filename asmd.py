#!/usr/bin/env python

"""asmd.py: ACSPO SPT Mask Differences
   -----------------------------------

   This program helps visualize questionable areas within the SPT mask
   Questionable areas are defined as pixels where the ACSPO mask is labeled
   non-ocean while the SPT mask has labeled the pixels ocean.  These areas
   are visualized by running a connected component algorithm on the questionable
   areas and created images centered around the areas.

   overlay description:
   overlay = 0: (acspo_mask==clear) & (spt_mask==clear)
   overlay = 1: (acspo_mask==clear) & (spt_mask==cloud)
   overlay = 2: (acspo_mask==clear) & (spt_mask==front)
   overlay = 3: (acspo_mask==cloud) & (spt_mask==clear) 
   overlay = 4: (acspo_mask==cloud) & (spt_mask==cloud)
   overlay = 5: (acspo_mask==cloud) & (spt_mask==front)"""


import os
import sys

import netCDF4
from skimage import measure
from scipy import ndimage
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as cl
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
        self.date_time = path.split('/')[-1][22:37]
        

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

    def generate_save_folder(self,loc):
        save_folder = os.path.join(loc,self.date_time)

        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        self.save_folder = save_folder

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

        self.height, self.width = layers['albedo'].shape
    	return layers


def generate_saveloc(up, down, left, right, save_folder):
    crop_name = "%d_%d_%d_%d" % (up, down, left, right)
    print crop_name
    return os.path.join(save_folder,crop_name)


if len(sys.argv) < 3:
	print "usage: python asmd.py <granule_path> <save_path>"
	sys.exit()

granule_path = sys.argv[1]
save_path = sys.argv[2]

granule = Granule(granule_path)
granule.generate_save_folder(save_path)
layers = granule.get_layers()

#create custom colorbar for overlay
colors = ['#000000', '#4DAF4A', '#999999', '#E41A1C', '#377EB8','#FF7F00']
cmap_overlay = cl.ListedColormap(colors)

#compute connected components
labels, n_lbls = measure.label(layers['questionable'],background=0, return_num=True)

#compute centers of clusters and left right up and down of cluster and save to file
centers = {}
inds_left = {}
inds_right = {}
inds_down = {}
inds_up = {}

padding = 50
for lbl in range(n_lbls):
    mask = labels == lbl
    centers[lbl] = ndimage.measurements.center_of_mass(labels,mask)
    inds_row, inds_col =  np.where(mask)
    
    inds_left[lbl] = np.min(inds_col)
    inds_right[lbl] =np.max(inds_col)
    inds_up[lbl] = np.min(inds_row)
    inds_down[lbl] = np.max(inds_row)

    
    # add padding, check boundaries
    if(inds_left[lbl] >= padding):
        inds_left[lbl] -= padding
    else:
        inds_left[lbl] = 0 

    if(inds_right[lbl] < granule.width - padding):
        inds_right[lbl] += padding
    else:
        inds_right[lbl] = granule.width - 1 

    if(inds_up[lbl] >= padding):
        inds_up[lbl] -= padding
    else:
        inds_up[lbl] = 0 

    if(inds_down[lbl] < granule.height - padding):
        inds_down[lbl] += padding
    else:
        inds_down[lbl] = granule.height - 1 
    
    crop = np.s_[inds_up[lbl]:inds_down[lbl], inds_left[lbl]:inds_right[lbl]]
    saveloc = generate_saveloc(inds_up[lbl],inds_down[lbl], inds_left[lbl],inds_right[lbl], granule.save_folder)
    
    fig = plt.figure()
    fig.suptitle("Diagnostic at Crop\n%d : %d, %d : %d" % (inds_up[lbl],inds_down[lbl], inds_left[lbl],inds_right[lbl]))

    ax1 = plt.subplot(221)
    ax1.set_title("Sea Surface Temperature")
    ax1.axis('off')
    img1 = ax1.imshow(layers['sst'][crop],vmin=271.15,vmax=307)
    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    ax2 = plt.subplot(222,sharex=ax1,sharey=ax1)
    ax2.set_title('Albedo')
    ax2.axis('off')
    img1 = ax2.imshow(layers['albedo'][crop], vmin=0, vmax=10, cmap='gray')
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    ax2 = plt.subplot(223,sharex=ax1,sharey=ax1)
    ax2.set_title('SST Gradient Magnitude')
    ax2.axis('off')
    img1 = ax2.imshow(layers['gradient'][crop],vmin=0, vmax=0.5, cmap='gray')
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    norm = cl.BoundaryNorm(np.arange(0,6,1),cmap_overlay.N )
    ticklabels = ['Agreed Clear' ,'SPT Cloud','Clear Front','SPT Clear', 'Agreed Cloud','Cloud Front']
    ax2 = plt.subplot(224,sharex=ax1,sharey=ax1)
    ax2.set_title('Questionable ')
    ax2.axis('off')
    img1 = ax2.imshow(layers['overlay'][crop],cmap=cmap_overlay,vmin=-0.5,vmax=5.5)
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1,ticks=np.linspace(0,5,6))
    cbar1.set_ticklabels(ticklabels)

    plt.savefig(saveloc, bbox_inches='tight')
    plt.close()
