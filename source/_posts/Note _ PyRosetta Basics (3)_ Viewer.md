title: 'PyRosetta Basics (3): Viewer'
date: '2025-08-21 12:43:37'
updated: '2025-08-21 12:43:40'
tags:
  - PyRosetta
  - protein design
  - biochemistry
  - note
  - protein
categories:
  - protein design
---
# PyRosetta Basics (3): Viewer
## Introduction
The pyrosetta.distributed Viewer quickly renders `.pdb` files, dynamically instantiating Pose objects if required for certain visualization modules (`viewer.set*`). 

So when adding visualization modules to the Viewer or using presets, passing Pose or PackedPose objects to the Viewer is suggested for quicker rendering. 

If a Pose object or list, tuple, or set of Pose objects are provided to the Viewer, the Viewer will dynamically update upon Pose conformational changes by calling the view.show() method or equivalently view(). 

The Viewer applies visualization modules in the same order they are added (from left to right), so layering different styles (and ResidueSelectors) on top of one another becomes possible. 

<font style="color:#DF2A3F;">WARNING</font>: `pyrosetta.distributed.viewer`runs in a jupyter notebook.

## Initialization & load
```python
!pip install pyrosettacolabsetup
import pyrosettacolabsetup; pyrosettacolabsetup.install_pyrosetta()
import pyrosetta; 
pyrosetta.init()

import glob
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import pyrosetta.distributed.viewer as viewer
import sys
```

```python
my_pose = pyrosetta.io.pose_from_file("inputs/3EK4.pdb")
```

## Configuring `ditributed`for visualizing ligands and non-canonical residues
The user must have already initialized PyRosetta providing `.params` for any ligands and non-canonical residues in the input `Pose`objects /`PackedPose`objects/`.pdb` files.

Otherwise `pyrosetta.distributed` automatically initializes PyRosetta with default command line options.

E.g.:

```python
flags = """
-auto_setup_metals 1
-detect_disulf 1
"""
# Display metal: True
# Display disulfide bonds: True
pyrosetta.distributed.init(flags)
```

## Instantiation of `distributed.viewer`
```python
view = viewer.init(pose, window_size=(800, 600))
view()
```

## `distributed.viewer`methods
The pyrosetta.distributed Viewer quickly renders `.pdb` files, dynamically instantiating Pose objects if required for certain visualization modules (`viewer.set*`). 

```python
# All Available `viewer` objects:
viewer.__all__

'''
>>>
['expand_notebook',
 'init',
 'presets',
 'setBackgroundColor',
 'setDisulfides',
 'setHydrogenBonds',
 'setHydrogens',
 'setStyle',
 'setSurface',
 'setZoom',
 'setZoomTo']
 
'''
```

1. visualization modules: `viewer.set*`
2. visualization presets: `viewer.presets`

```python
viewer.presets.__all__
# >>> ['coreBoundarySurface', 'ligandsAndMetals']
```

3. `pyrosetta.distributed.viewer.expand_notebook()` expands the Jupyter notebook cell width to fit your internet browser

