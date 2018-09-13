Crystal-proj
============


Crystal-proj is a python script that draws the projection of a crystal structure on a given plane and the dichromatic pattern.

## Requirements
* Python 2.7 with Matplotlib, Numpy, Tkinter, PIL (available through [Enthought Canopy](https://store.enthought.com/downloads/) or [Anaconda](http://continuum.io/downloads) distributions for instance).
* Mac users need to install [Active Tcl](http://www.activestate.com/activetcl/downloads).
* Run the script python /my-path-to-the-script/crystal-proj.py

## User guide

### Interface
![img1](/img1.png?raw=true)

### Procedure
* Enter the crystal structure or import it from the menu bar. The structure can be modified/added by modifying the structure.txt file. The format is: name a b c alpha beta gamma space-group. 
* Enter the indices of the plane to draw
* Enter the layers number (i.e. the number of parallel planes)
* Enter the size of the projectioned plane
* Click calculate. You will be prompted to select the appropriate crystal structure. Text files such as cfc.txt describe the position of atoms in the unit cell. The format is: 

atom1 x1 y1 z1

...

atomn xn yn zn

v1x v1y v1z

...

vkx vky vkz

with atomn is nth atom in the cell at coordinate xn, yn, zn and vkx,vky,vkz is the kth translation vector defining all the positions of atoms in the cell.
* Click on the draw button.

* Tick "atoms" to see the different atoms (label with different numbers)

* Tick "label" to see the coordinates

![img2](/img2.png?raw=true)

* Change the marker size with the scale button and shape using "circle/square" button

* Zoom in, move and save the figure using the bottom tool bar

* Press 3D button to draw a 3D draggable/rotating view of the crystal.

![img3](/img3.png?raw=true)

* To draw the corresponding dichromatic pattern check the "DSC" button, recalculate (press "Calculate") and then "Draw". The dichromatic pattern is drawn according the usual convention: black and white crystals correspond to two different grains and different symbols to different layers (note that you can plot up to 7 different layers with different symbols)

![img4](/img4.png?raw=true)

* Tick the "Coincidence sites" button to draw coincident site positions. The precision by which is defined coincidence  can be tuned. It is expressed in % of the largest lattice parameter.

![img5](/img5.png?raw=true)

* Tick the "Layers" button to show the layer number on the dichromatic pattern. From layers '0' to '5', markers are circles, squares, triangles, stars and hexagons, respectively. After they are all circles. 
