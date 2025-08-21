title: 'PyRosetta Basics (2): PyMOL Mover'
date: '2025-08-21 12:41:34'
updated: '2025-08-21 12:41:37'
tags:
  - PyRosetta
  - biochemistry
  - protein design
  - protein
  - note
categories:
  - protein design
---
# PyMOL Mover
```python
from pyrosetta import *
init()
#import pyrosetta
#pyrosetta.init()
```

```python
my_pose = pose_from_pdb("inputs/5tj3.pdb")
```

## Resource
[https://www.pyrosetta.org/documentation/pymol_mover-tutorial](https://www.pyrosetta.org/documentation/pymol_mover-tutorial)

## `PyMOLMover` class
The PyMOLMover class will let us send information from PyRosetta to PyMOL for quick visualization. 

### Setup PyMOL
#### Get Configuration file: `<font style="color:rgb(33, 33, 33);background-color:rgba(0, 0, 0, 0.06);">PyMOLPyRosettaServer.py</font>`(in Linux)
<font style="color:rgb(33, 33, 33);background-color:rgba(0, 0, 0, 0.06);">PyMOLPyRosettaServer.py</font><font style="color:rgb(33, 33, 33);"> is found in the main directory of PyRosetta.</font>

<font style="color:rgb(33, 33, 33);">The path to PyRosetta can be discovered in the init information.</font>

<font style="color:rgb(33, 33, 33);">For example:</font>

```bash
Requirement already satisfied: pyrosettacolabsetup in /home/michael2003/anaconda3/envs/DESIGN_PROTEIN/lib/python3.9/site-packages (1.0.9)
┌──────────────────────────────────────────────────────────────────────────────┐
│                                 PyRosetta-4                                  │
│              Created in JHU by Sergey Lyskov and PyRosetta Team              │
│              (C) Copyright Rosetta Commons Member Institutions               │
│                                                                              │
│ NOTE: USE OF PyRosetta FOR COMMERCIAL PURPOSES REQUIRE PURCHASE OF A LICENSE │
│         See LICENSE.PyRosetta.md or email license@uw.edu for details         │
└──────────────────────────────────────────────────────────────────────────────┘
PyRosetta-4 2025 [Rosetta PyRosetta4.conda.ubuntu.cxx11thread.serialization.Ubuntu.python39.Release 2025.24+release.8e1e5e54f047b0833dcf760a5cd5d3ce94d63938 2025-06-06T09:20:57] retrieved from: http://www.pyrosetta.org
core.init: Checking for fconfig files in pwd and ./rosetta/flags
core.init: Rosetta version: PyRosetta4.conda.ubuntu.cxx11thread.serialization.Ubuntu.python39.Release r403 2025.24+release.8e1e5e54f0 8e1e5e54f047b0833dcf760a5cd5d3ce94d63938 http://www.pyrosetta.org 2025-06-06T09:20:57
core.init: Rosetta extras: [cxx11thread, serialization]
core.init: command: PyRosetta -ex1 -ex2aro -database /home/michael2003/anaconda3/envs/DESIGN_PROTEIN/lib/python3.9/site-packages/pyrosetta/database
basic.random.init_random_generator: 'RNG device' seed mode, using '/dev/urandom', seed=-481462242 seed_offset=0 real_seed=-481462242 thread_index=0
basic.random.init_random_generator: RandomGenerator:init: Normal mode, seed=-481462242 RG_type=mt19937
```

In this example, path to PyRosetta is:

`/home/michael2003/anaconda3/envs/DESIGN_PROTEIN/lib/python3.9/site-packages/pyrosetta/`

#### Copy & edit the configuration file (in Windows)
I'm using PyMol on Windows, and pyRosetta on WSL Ubuntu.

Therefore, I have to edit the <font style="color:rgb(33, 33, 33);background-color:rgba(0, 0, 0, 0.06);">PyMOLPyRosettaServer.py </font>file to appoint a correct IP address.

[https://learn.microsoft.com/en-us/windows/wsl/networking](https://learn.microsoft.com/en-us/windows/wsl/networking)

> <font style="color:rgb(22, 22, 22);">A program running inside a Linux distribution via WSL2 (instance) wants to know the Windows host's IP address, so that a Linux program can connect to a Windows host server program.</font>
>

```bash
ip route show | grep -i default | awk '{ print $3}'
```

And a typical output might look like:

```powershell
>>> 172.24.112.1
```

Copy the <font style="color:rgb(33, 33, 33);background-color:rgba(0, 0, 0, 0.06);">PyMOLPyRosettaServer.py </font><font style="color:rgb(33, 33, 33);">  from Linux to Windows (to D:).</font>

```bash
cp <Path-to-PyRosetta>/PyMOLPyRosettaServer.py /mnt/d/
```

And edit the file according to the comments within it. (The IP configuration is in the last few lines.) Appoint the right remote IP, and appoint an available port (e.g. `65000`).

#### Apply the <font style="color:rgb(33, 33, 33);background-color:rgba(0, 0, 0, 0.06);">PyMOLPyRosettaServer.py </font><font style="color:rgb(33, 33, 33);"> to PyMOL (in Windows)</font>
a. Drag the file to PyMOL window.

b. Or, run the following commands in the PyMOL command line.

```powershell
cd D: #Because I save PyMOLRosettaServer.py there.
run PyMOLRosettaServer.py
```

c. Or, save the commands above to `pymolrc.pml`; thus it runs automatically each time pyMOL launches  
(The `pymolrc.pml`file can be found from GUI menu: `File` >> `Edit pymolrc`)

### Instantiation of `PyMOLMover`
```python
pmm = PyMOLMover("172.24.112.1", 65000)
```



### `PyMOLMover`methods
#### Keep history
The method `keep_history`, if set to True, allows you to load in structures with the same name into states of the same object in PyMOL. 

This is the starting point for creating a PyMOL movie of your structure.

```python
pmm.keep_history(True) # This keeps history of all poses sent to PyMOL
```

#### Send a `pose` object to PyMOL
```python
pmm.apply(my_pose)
# After running this, my_pose is expected to show up in the PyMOL window
```

See `AddPyMOLObserver`class for information about automatic updating the _pose_ during _move._

#### Send energy (energy across every residue) of a `pose`object to PyMOL
```python
# assume that `pmm` is an instance of `PyMOLMover`
# assume that `my_pose` is an instance of `Pose`
pmm.send_energy(my_pose) # visualize the distribution of energy across the structure
pmm.send_energy(my_pose, fa_sol) # visualize the distribution of solvation energy (of full atom energy landscape) across the structure
```

#### Send properties of the `pose`object to PyMOL
```python
# if you have scored the pose first, you can:
pmm.send_hbonds(my_pose) #show all hydrogen bonds in PyMOL
```

## `AddPyMOLObserver`class
<font style="color:rgba(0, 0, 0, 0.87);">The observer is configured to execute a </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">PyMOLMover.apply()</font>`<font style="color:rgba(0, 0, 0, 0.87);"> every time a change is observed in the pose coordinates. The </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">True</font>`<font style="color:rgba(0, 0, 0, 0.87);"> is a flag to ensure that PyMOL keeps a history of the moves.</font>

<font style="color:rgba(0, 0, 0, 0.87);">The observer should be added before excuting the </font>_<font style="color:rgba(0, 0, 0, 0.87);">move(s).</font>_

```python
observer = pyrosetta.rosetta.protocols.moves.AddPyMOLObserver(my_pose, True)
my_mover.apply(my_pose)
```

