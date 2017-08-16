#ACSPO SPT Mask Differences

This program helps visualize questionable areas within the SPT mask
Questionable areas are defined as pixels where the ACSPO mask is labeled
non-ocean while the SPT mask has labeled the pixels ocean.  These areas
are visualized by running a connected component algorithm on the questionable
areas and created images centered around the areas.

The necessary libraries needed to run this program can be found in the file [requirements.txt](requirements.txt)
Note the program has not been tested with any other versions of the provided libraries.

To run the program:

`python asmd.py path-to-granule path-save-folder`