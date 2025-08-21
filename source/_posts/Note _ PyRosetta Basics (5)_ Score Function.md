# PyRosetta Basics (5): Score Function
```python
!pip install pyrosettacolabsetup
import pyrosettacolabsetup; pyrosettacolabsetup.install_pyrosetta()
import pyrosetta; pyrosetta.init()
import pyrosetta
pyrosetta.init()
```

## Introduction
Rosetta scores the energy of a `pose`object with a score function.

The lower the enegy, the more stable the pose.

The score function is an object of `pyrosetta.ScoreFunction` class.

The total score is a sum of per-residue score.

The unweighted value of energy terms (of each residue) is stored in the `pose`object, and can be visited via `Pose.energies`method. (See [Note | PyRosetta Basics (1): Pose](https://www.yuque.com/yuqueyonghu2r84jv/zmy8pd/ktcfk65y26wfqciw))

The score function stores a series of weights correspond to each energy term, and calculate the total score.

[Energy Terms In Rosetta](https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types)

[Centroid score terms and score functions in Rosetta](https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/centroid-score-terms)

[Scoring PyRosetta 4.0 document](https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.scoring.html)

## Instantiation of `ScoreFunction`
[Centroid score terms and score functions in Rosetta](https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/centroid-score-terms)

### Empty score function
```python
sfxn = ScoreFunction() #sfxn is an empty `score_function` object

# and we can manually set the weights
sfxn.set_weight(fa_atr, 1.0) 
sfxn.set_weight(fa_rep, 1.0)
```

### `teaching.get_score_function`method
```python
from pyrosetta.teaching import *
sfxn = get_score_function(True) 
# get_score_function(is_fullatom: bool) in pyrosetta.teaching namespace
# True: default ref2015 all-atom energy function
# False: default centroid score function
```

```python
cen_sfxn = pyrosetta.create_score_function("score0") 
# a centroid score function
# score function used in the first stage of the ClassicAbInitio protocol
# https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/centroid-score-terms
```

### `<font style="color:rgb(33, 33, 33);">get_fa_scorefxn</font>`method (the default ref2015 function)
<font style="color:rgb(33, 33, 33);">The </font>`<font style="color:rgb(33, 33, 33);">get_fa_scorefxn()</font>`<font style="color:rgb(33, 33, 33);"> method is a wrapper of a method within </font>`<font style="color:rgb(33, 33, 33);">ScoreFunctionFactory</font>`<font style="color:rgb(33, 33, 33);"> which returns the current standard </font>`<font style="color:rgb(33, 33, 33);">ScoreFunction</font>`<font style="color:rgb(33, 33, 33);">.</font>

```python
scorefxn = get_fa_scorefxn()
scorefxn(my_pose)
```

### `create_score_function`(specify a `weights_tag`)
```python
cart_sf = create_score_function("ref2015_cart") # Cartesian for Regional FastRelax
```

## Score a pose
```python
my_pose = pyrosetta.pose_from_pdb("inputs/6Q21_A.pdb") # instantiate a `pose` object
my_score = sfxn(my_pose)
# my_score is the total score; type is float
```

## `ScoreFunction` methods
### Show energy terms in detail
```python
my_sfxn.show(my_pose)
'''
>>>
core.scoring: 
------------------------------------------------------------
 Scores                       Weight   Raw Score Wghtd.Score
------------------------------------------------------------
 fa_atr                       1.000   -1039.246   -1039.246
 fa_rep                       0.550    1193.837     656.611
 fa_sol                       1.000     682.582     682.582
 fa_intra_rep                 0.005     700.419       3.502
 fa_intra_sol_xover4          1.000      46.564      46.564
 lk_ball_wtd                  1.000     -14.597     -14.597
 fa_elec                      1.000    -195.387    -195.387
 pro_close                    1.250      97.210     121.513
 hbond_sr_bb                  1.000     -41.656     -41.656
 hbond_lr_bb                  1.000     -28.352     -28.352
 hbond_bb_sc                  1.000     -13.111     -13.111
 hbond_sc                     1.000      -7.771      -7.771
 dslf_fa13                    1.250       0.000       0.000
 omega                        0.400      41.525      16.610
 fa_dun                       0.700    1296.642     907.650
 p_aa_pp                      0.600     -25.496     -15.298
 yhh_planarity                0.625       0.000       0.000
 ref                          1.000      47.114      47.114
 rama_prepro                  0.450     197.781      89.002
---------------------------------------------------
 Total weighted score:                     1215.729
 '''
```

### Show Intrinsic weights of a scoring function
```python
scorefxn = get_fa_scorefxn() 
# First create object by calling wrapped methods in `get_fa_scorefxn()`
scorefxn.weights()[fa_atr]
# Then get data from the `scorefxn` object
```

Or, use `ScoreFunction`class:

```python
scorefxn = ScoreFunction()
# First create object by calling wrapped methods in `get_fa_scorefxn()`
scorefxn.weights()[fa_atr]
# Then get data from the `scorefxn` object
```

<font style="color:#DF2A3F;">WARNING</font>: The following usage is incorrect:

```python
get_fa_scorefxn().weights()[fa_atr] #incorrect; returns 0.0
```

Although the following lines do not cause mistakes, chain calling is a depreciated practice.

```python
ScoreFunction().weights()[fa_atr]
# This creates an object `ScoreFunction`, 
#  and then shows the default weight of the term `fa_atr`
# Works fine, but depreciated.
```

## `EMapVector`class
### Introduction
`EMapVector`class is an auxillary class for energy.

An instance of the`EMapVector`class is a vector that can be used to stores energies.

### Example: Analyzing energy between residues
```python
my_pose = pyrosetta.toolbox.pose_from_rcsb("1YY9") 
# instantiate a `pose` object

residue_102_id = my_pose.pdb_info().pdb2pose("D", 102) 
# get pose id of residue #102(PDB) on chain D
residue_408_id = my_pose.pdb_info().pdb2pose("A", 408) 
# get pose id of residue #408(PDB) on chain A

residue_102 = my_pose.residue(residue_102_id)
residue_408 = my_pose.residue(residue_408_id)
# get `residue` objects

my_emap = EMapVector() 
# Instantiate an empty `e_map_vector` object

sfxn = get_score_function(True) 
# Create a score function `sfxn` using
#  pyrosetta.teaching.get_score_function(); ref2015 scoring function

sfxn.eval_ci_2b(residue_102, residue_408, my_pose, my_emap)
# Analyzing energy between residues #102(chain D, PDB) and #408(chain A, PDB)
# The outputs are stored in my_emap

print(emap[fa_atr])
print(emap[fa_rep])
print(emap[fa_sol])
# visit and print the energy terms that are stored in emap
```



## `SwitchResidueTypeSetMover`class
An instance of the`SwitchResidueTypeSetMover`class is a `mover`that switches the energy landscape of a pose.

+ "centroid": a simplified, smooth energy landscape; easy to navigate; fast
+ "fa_standard":  a full-atom landscape; rough and sharp; accurate

PyRosetta uses different energy terms and different scoring functions for these two energy landscape.

```python
# assume that we have a `pose` object my_pose
c_switch = SwitchResidueTypeSetMover("centroid") 
# `c_switch` is an object that can switch a pose to low-resolution (centroid) energy landscape
switch.apply(my_pose) # apply switch to my_pose
```

```python
fa_switch = SwitchResidueTypeSetMover("fa_standard")
# `fa_switch` is an object that can switch a pose to high-resolution (full atom) energy landscape
switch.apply(my_pose)
```

## Notes on Energy Terms `ScoreType`
### Introduction
namespace: `pyrosetta.rosetta.core.scoring`.

E.g. 

```python
scorefxn = get_fa_scorefxn()
scorefxn(pose)
energies = pose.energies()
#print(energies.residue_total_energies(49))
print(energies.residue_total_energies(49)[pyrosetta.rosetta.core.scoring.fa_dun])
```

### Notes on Abbreviations
`fa_*`indicates that this term is measured in a "full atom" energy landscape.

`bb`means backbone.

`sc`means side chain.

`vdw`for Van der Waal.

`fa_dun`for Dunbrak Energies (of the side chain)



