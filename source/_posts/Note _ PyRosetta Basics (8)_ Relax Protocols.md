title: 'PyRosetta Basics (8): Relax Protocols'
date: '2025-08-21 12:39:12'
updated: '2025-08-21 12:39:17'
tags:
  - protein
  - biochemistry
  - protein design
  - PyRosetta
  - note
categories:
  - protein design
---
# PyRosetta Basics (8): Relax Protocols
## Initialization(Same as the "Packing" chapter)
```python
#Python
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *

#Core Includes
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.simple_metrics import metrics
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.core.select.movemap import *

#Protocol Includes
from pyrosetta.rosetta.protocols import minimization_packing as pack_min
from pyrosetta.rosetta.protocols import relax as rel
from pyrosetta.rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.loops import *
from pyrosetta.rosetta.protocols.relax import FastRelax
```

```python
init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false -relax:default_repeats 2 -no_fconfig')
```

## Introduction to Relax
### Overview
> Relax is the main protocol for simple all-atom refinement of structures in the Rosetta force-field. Relax does not do extensive refinement and only searches the local conformational space around the starting structure. Relax is thus often used in conjunction with more aggressive sampling protocols like fragment assembly (abinitio) and loop modelling. To evaluate different conformations based on their Rosetta all-atom score one usually has to apply relax.
>
> It can also read centroid models, in which case it will convert the model into a fullatom model and pack the sidechains. Relax does not carry out any extensive refinement and only searches the local conformational space neighbourhood.
>
> It is further advisable to apply relax only to previously idealized structures. Idealization avoids that score differences arise due to non-ideal geometry (e.g., at the position of former chain-breaks introduced during an aggressive sampling stage and removed by loop closing).
>
> [https://docs.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax](https://docs.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax)
>

### `ClassicRelax`
Depreciated. Small move + shear move.

### `FastRelax`(default)
5 cycles of {[packing + minimization] * n_ramping}.

> Packing: optimizing side chain conformations without editing side chains.
>
> Minimization: find local energy minima.
>
> n_ramping: steps of ramping up `fa_rep` (in one single cycle).
>

In each of the 5 cycles, the repulsion energy (`fa_rep`) are firstly set to a small number, and then gradually increase to a realistic number in n-step ramping.

Pseudo code:

```python
# Assume that we have a `Pose` object called `my_pose`
for i in range(0,5): # "5" by default, can be modified
    fa_rep = 0.1 * fa_rep
    for j in range(1, n_ramp):
        fa_rep = Increase(fa_rep)
        Packing(my_pose)
        Minimizing(my_pose)
```

### `CentroidRelax`(for massive screening)
1. It uses centroid score functions
2. It ramps up various energy terms while minimizing the pose
3. It is rough but fast, thus is suitable for high-throughput structure design   
(as a pre-processing & evaluation for further `FastRelax`).





## `FastRelax` class
### Instantiation & Configurations
```python
my_fr = FastRelax() # Instantiation a `FastRelax` object named `my_fr`

my_scorefxn = get_score_function() # A default ref2015 full atom score function

my_fr.set_scorefxn(my_scorefxn)

#FastRelax takes a very long time, 
#but we can decrease the amount of minimization cycles we use:
#(Only recommended for cartesian)
# fr.max_iter(100)
```

### Apply `FastRelax`to a `Pose`object
```python
# Assume that we have a `Pose` object named `my_pose`
my_fr.apply(my_pose)
```

## Regional Relax: `MoveMapFactory`class
### Introduction to `MoveMapFactory` class
`MoveMapFactory`is the prototype of `MoveMap`.

`MoveMap` instructs a mover. 

The default is to have everything OFF first, and turn specific things on.



Regional Relax requires:

1. pre-defined residue selectors
2. the pose

So residue selectors are designed to be passed to the `MoveMapFactory`.

[Note | PyRosetta Basics (6): Movers for Backbone Folding](https://www.yuque.com/yuqueyonghu2r84jv/zmy8pd/gfkfxof6pygtuwkr)

[MoveMapFactories (RosettaScripts)](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/MoveMapFactories/MoveMapFactories-RosettaScripts)

### Instantiation & configuration
```python
# Establish Residue Selector
cdr_selector = CDRResidueSelector()
cdr_selector.set_cdr(h1)
# h1 is a Enum imported when we imported the antibody namespace

# Instantiate a `MoveMapFactory` object named "my_mmf"
my_mmf = MoveMapFactory()
# Setup `my_mmf`
my_mmf.add_bb_action(mm_enable, cdr_selector)
my_mmf.add_chi_action(mm_enable, cdr_selector) 
# mm_enable and mm_disable are Enums (numbered variables) that 
#  come when we import the MMF.
```

### Check the `MoveMapFactory` settings
Based on the `MoveMapFactory`, we can create a `MoveMap`object from a `Pose`object, and thus we can check the settings.

```python
# Instantiate a `MoveMap` object from `my_mmf`
# We need to pass in a `Pose` object (`my_pose`)
my_mm  = my_mmf.create_movemap_from_pose(my_pose)

# Inspect my_mm
print(my_mm)
```

### Basic Setup of Regional `FastRelax`
```python
# Assume that we have a `FastRelax` object called `my_fr`
# Assume that we have a `MoveMapFactory` object called `my_mmf`
my_fr.set_movemap_factory(my_mmf)

# Assume that we have a `TaskFactory` object called `pack_cdrs_and_neighbors_tf`
my_fr.set_task_factory(pack_cdrs_and_neighbors_tf)

# Assume that we have a default ref2015 full atom score function
my_fr.set_scorefxn(my_scorefxn)

# Assume that we have a `Pose` object called `my_pose`
# Finally, we can apply the `my_fr` to `my_pose`:
my_fr.apply(my_pose)
```

### Optimization
Basic Setup may cause large conformation shifts that we don't want.

We can tackle this problem by using a different score function, OR by customizing the `FoldTree`.

(Not both)

#### Cartesian-space refinement `ref2015_cart`
```python
# Create score function specified for Cartesian
cart_sf = create_score_function("ref2015_cart")

# Assume that we have a `MoveMapFactory` object called `my_mmf`
my_mmf.set_cartesian(True) # Turn "Cartesian mode" on.

# Assume that we have a `FastRelax` object called `my_fr`
my_fr.set_movemap_factory(my_mmf)

# Use `cart_sf` as the score function for FastRelax
my_fr.set_scorefxn(cart_sf)
my_fr.cartesian(True)

#This is a general recommendation for cartesian minimization - it lowers the number of maximum cycles.
# More than this only increases time of protocol, but has little effect on energies/structure
fr.max_iter(200)
```

These settings are Cartesian-specifc. They are turned off by default. Turn them off under other circumstances.

```python
my_mmf.set_cartesian(False)
my_fr.cartesian(False)
my_fr.max_iter(0) #Reset to default 
```

#### Customizing `FoldTree` (A Classic way)
<font style="color:#DF2A3F;">WARNING</font>: If Cartesian has been turned on, we should turn Cartesian off before using this solution.

```python
# Assume that we have a `Pose` object named `my_pose`, and it is an antibody
ab_info = AntibodyInfo(my_pose) # Get antibody-specific information
my_ft = FoldTree() # Instantiate an empty FolTree

# Get CDR H1 loop parameters
# >>> #
start = ab_info.get_CDR_start(h1, my_pose)
stop =  ab_info.get_CDR_end(h1, my_pose)
cutpoint = int((stop-start)/2) + start
# <<< #

# Define the loop
cdr_loop = Loop(start, stop, cutpoint)
cdr_loops = Loops()
cdr_loops.add_loop(cdr_loop)

# Setup of `my_ft`
fold_tree_from_loops(my_pose, cdr_loops, my_ft) 

my_pose.fold_tree(my_ft)
original_ft = my_pose.fold_tree() # Backup the original fold tree
add_cutpoint_variants(my_pose)

# Create a default full atom score function
scorefxn = get_score_function()

# Add chainbreak term to `scorefxn` so we don't get wacky stuff.  
# This term helps keep the peptide closed during bb movement.
scorefxn_ch = scorefxn # Copy the default score function
scorefxn_ch.set_weight(rosetta.core.scoring.chainbreak, 100) # Add chainbreak term

# Assume that we have a `FastRelax` object called `my_fr`
my_fr.set_scorefxn(scorefxn_ch)

# Assume that we have a `MoveMap` object (with basic configurations) called `my_mmf`
my_fr.set_movemap_factory(my_mmf)

my_fr.max_iter(0) #Reset to default 
# if it's 0, then we don't set it in the MinMover that FastRelax runs

# Start FastRelax
my_fr.apply(my_pose)

#Reapply the original fold tree
my_pose.fold_tree(original_ft)
```

