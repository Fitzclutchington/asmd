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
        path: a string of the path to the granule file
        date_time = a string of the date and time of the granule
        day = an boolean denoting whether the granule was taken at day or night
        height = an integer representing the number of rows of the granule
        width = an integer representing the number of columsn in the granule
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

    def	_get_acspo_values(self, acspo_flags):
    	mask_bits = 192
        day_bit = 2
        land_bit = 4

    	mask = np.bitwise_and(acspo_flags,mask_bits).astype(bool)
        land_mask = np.bitwise_and(acspo_flags,land_bit).astype(bool)
        self.day = np.bitwise_and(acspo_flags,day_bit).astype(bool)[0,0]

    	return (mask, land_mask)

    def _compute_overlay(self, acspo_mask, spt_mask):
        overlay = np.zeros(acspo_mask.shape).astype(np.uint8)
        overlay[(acspo_mask==0) & (spt_mask==0)] = 0 
        overlay[(acspo_mask==0) & (spt_mask==4)] = 1
        overlay[(acspo_mask==1) & (spt_mask==0)] = 2 
        overlay[(acspo_mask==1) & (spt_mask==3)] = 3
        overlay[(acspo_mask==1) & (spt_mask==4)] = 4 
    	return overlay

    def _compute_gradient(self, sst):
    	dX,dY = np.gradient(sst)
    	return np.sqrt(dX**2 + dY**2)

    def generate_save_folder(self,loc):
        save_folder = os.path.join(loc,self.date_time)

        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        self.save_folder = save_folder

    def _generate_land_layer(self,land_mask):
        land = np.zeros((self.height,self.width,4))
        r = 146/256.0
        g = 98/256.0
        b = 57/256.0
        land[land_mask] = [r,g,b,1]
        return land

    def get_layers(self):
    	"""open the cdf file and get the necessary variable
    	   return a dictionary with the layer names as keys and 
    	   matrices as value"""

    	layers = {}
        acspo_flags = self._read_var('acspo_mask')
        self.height, self.width = acspo_flags.shape
        acspo_mask, land_mask =  self._get_acspo_values(acspo_flags)

        layers['land'] = self._generate_land_layer(land_mask)

    	layers['sst'] = self._read_var('sst_regression')
    	
        spt_mask = self._read_var('spt_mask')
    	layers['spt'] = spt_mask  	
    	layers['overlay'] = self._compute_overlay(acspo_mask,spt_mask)
    	layers['gradient'] = self._compute_gradient(layers['sst'])
        layers['questionable'] = ((layers['overlay']==2) | (layers['overlay']== 4))

        if self.day:
            layers['albedo'] = self._read_var('albedo_chM7')
        else:
            reference = self._read_var('sst_reynolds')
            layers['ref_diff'] = layers['sst'] - reference

        
    	return layers


def generate_saveloc(up, down, left, right, save_folder):
    crop_name = "%d_%d_%d_%d" % (up, down, left, right)
    print crop_name
    return os.path.join(save_folder,crop_name)

def add_offset(ind, offset ,  alt, max_val):
    if not max_val:
        if(ind >= offset):
            ind -= offset
        else:
            ind = alt
    else:
        if(ind < (max_val - offset)):
            ind += offset
        else:
            ind = alt
    return ind

def generate_crop(mask):
    min_height = 400
    min_width = 400

    padding = 100
    inds_row, inds_col =  np.where(mask)
    
    inds_left = np.min(inds_col)
    inds_right =np.max(inds_col)
    inds_up = np.min(inds_row)
    inds_down = np.max(inds_row)

    height_crop = inds_down - inds_up
    width_crop = inds_right - inds_left

    # add padding, check boundaries
    if width_crop >= min_width:
        inds_left = add_offset(inds_left, padding, 0, 0)
        inds_right = add_offset(inds_right, padding, granule.width-1, granule.width)


    else:
    
        # add padding to the left
        width_remaining = min_width - width_crop
        inds_left = add_offset(inds_left, width_remaining/2, 0, 0)
   
        # add padding to the right
        width_remaining = min_width - (inds_right - inds_left)
        inds_right = add_offset(inds_right, width_remaining, granule.width-1, granule.width)
     
        # add remaining pdding to left if right was greater than width of matrix
        width_remaining = min_width - (inds_right - inds_left)
        inds_left = add_offset(inds_left, width_remaining, 0, 0)
      

    if height_crop >= min_height:
        inds_up = add_offset(inds_up, padding, 0, 0)
        inds_down = add_offset(inds_down, padding, granule.height-1, granule.height)


    else:

        # add padding to the top
        height_remaining = min_height - height_crop
        inds_up = add_offset(inds_up, height_remaining/2, 0, 0)
        
        # add padding to the bottom
        height_remaining = min_height - (inds_down - inds_up)
        inds_down = add_offset(inds_down, height_remaining, granule.height-1, granule.height)
        
        # add remaining pdding to top if bottom was greater than height of matrix
        height_remaining = min_height - (inds_down - inds_up)
        inds_up = add_offset(inds_up,height_remaining,0, 0)
    

    return (inds_up, inds_down, inds_left, inds_right)




if len(sys.argv) < 3:
	print "usage: python asmd.py <granule_path> <save_path>"
	sys.exit()

granule_path = sys.argv[1]
save_path = sys.argv[2]

granule = Granule(granule_path)
granule.generate_save_folder(save_path)
layers = granule.get_layers()

#create custom colorbar for overlay
colors = ['#000000', '#999999', '#E41A1C', '#377EB8','#FF7F00']
cmap_overlay = cl.ListedColormap(colors)

#compute connected components
labels, n_lbls = measure.label(layers['questionable'],background=0, return_num=True)

#compute centers of clusters and left right up and down of cluster and save to file
centers = {}


padding = 50
for lbl in range(n_lbls):
    mask = labels == lbl
    centers[lbl] = ndimage.measurements.center_of_mass(labels,mask)
    
    inds_up, inds_down, inds_left, inds_right = generate_crop(mask)

    crop = np.s_[inds_up:inds_down, inds_left:inds_right]
    saveloc = generate_saveloc(inds_up,inds_down, inds_left,inds_right, granule.save_folder)
    
    fig = plt.figure(figsize=(11,8))
    fig.suptitle("Diagnostic at Crop %d : %d, %d : %d" % (inds_up,inds_down, inds_left,inds_right))

    plt.ticklabel_format(useOffset=False)
    vmin = 271.15 if np.nanmin(layers['sst'][crop][layers['spt'][crop]!=3]) < 271.15 else np.nanmin(layers['sst'][crop][layers['spt'][crop]!=3])
    vmax = np.nanmax(layers['sst'][crop][layers['spt'][crop]!=3])

    if abs(vmin-vmax) < 1:
        vmin = vmin - 1
        vmax = vmax + 1

    ax1 = plt.subplot(221)
    ax1.set_title("Sea Surface Temperature", fontsize=10)
    ax1.axis('off')
    img1 = ax1.imshow(layers['sst'][crop],vmin=vmin,vmax=vmax)
    ax1.imshow(layers['land'][crop])
    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    ax2 = plt.subplot(222)
    ax2.axis('off')
    if granule.day:
        ax2.set_title('Albedo', fontsize=10)
        img1 = ax2.imshow(layers['albedo'][crop], vmin=0, vmax=10, cmap='gray')
    else:
        ax2.set_title('SST - Reference', fontsize=10)
        img1 = ax2.imshow(layers['ref_diff'][crop], vmin=-3, vmax=3)
    ax2.imshow(layers['land'][crop])
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    ax2 = plt.subplot(223)
    ax2.set_title('SST Gradient Magnitude', fontsize=10)
    ax2.axis('off')
    img1 = ax2.imshow(layers['gradient'][crop],vmin=0, vmax=0.5, cmap='gray')
    ax2.imshow(layers['land'][crop])
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1)

    ticklabels = ['Agreed Clear' ,'Clear Front','SPT Clear', 'Agreed Cloud','Cloud Front']
    ax2 = plt.subplot(224)
    ax2.set_title('Questionable', fontsize=10)
    ax2.axis('off')
    img1 = ax2.imshow(layers['overlay'][crop],cmap=cmap_overlay,vmin=-0.5,vmax=4.5)
    ax2.imshow(layers['land'][crop])
    div1 = make_axes_locatable(ax2)
    cax1 = div1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(img1, cax=cax1,ticks=np.linspace(0,4,5))
    cbar1.set_ticklabels(ticklabels)

    plt.savefig(saveloc)
    plt.close()
