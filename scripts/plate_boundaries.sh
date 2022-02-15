#!/bin/bash

#######################################################################################
# This script plots the coordinates of plate boundaries for one of the four databases: 
# (1) nuvel_plates.txt      : Nuvel Plate boundaries (Argus and Gordon, 1991)
# (2) gem_active_faults.txt : Global Earthquake Model (Pagani et al., 2018; Styronand Pagani, 2020)
# (3) PB2002_steps.txt      : Bird plate boundaries (Bird, 2003)
# (4) bird_gem_faults.txt   : Bird plate boundaries from the GEM database.

# The files (2), (3), and (4) can be generated using the notebook, make_world_builder_plate_boundaries.ipynb 
# by running the 2nd, 3rd, 4th cell (i.e., input the json file and use the faults to define each 'feature') 
# and the second-to-last cell (outputs the longitude, latitude in a text file such that each pair is in new line).
# Thereafter we remove the square brackets "[" and "]" in the output text files using vi and command "%s/\]//" 
# and "%s/\[//", and then use this bash file.
# To use the specific file names change the $file_name parameter for (2), (3), (4). 

# To plot nuvel plate boundaries, we use the file, nuvel_plates.txt, which already has line segments
# separated by ">" character.
########################################################################################

file_name=bird_gem_faults.txt 

# Plot coastlines on the globe with continents filled with gray color.
gmt pscoast -Rd -JN5i -Ggray -Swhite -Cwhite -Dc -K -P > plate_boundaries.ps

# Plot the gridlines and specify annotaions in the map.
gmt psbasemap -R -J -Ba60g60/a60g30 -P -O -K >> plate_boundaries.ps

# The following line is uncommented to plot the files output from the notebook. psxy plots the longitude,
# latitude points as filled black circles using psxy. 
# awk 'BEGIN {FS = ","}; {print $1, $2}' ${file_name} | gmt psxy -R -J -Sc1p -Wblack -P -O >> plate_boundaries.ps

# The following line plots the longitude, latitude segments separated by the character ">".
gmt psxy nuvel_plates.txt -R -J -W1p,black -m">" -P -O >> plate_boundaries.ps