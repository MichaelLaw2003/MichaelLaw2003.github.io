title: 'PyRosetta Basics (7): Movers for Side Chain Packing and Design'
tags:
  - PyRosetta
  - biochemistry
  - protein
  - note
  - protein design
categories:
  - protein design
---
# PyRosetta Basics (7): Movers for Side Chain Packing and Design
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

[https://docs.rosettacommons.org/docs/latest/full-options-list](https://docs.rosettacommons.org/docs/latest/full-options-list)

## Basic Concepts
### Terms of pose operations
1. "Folding": Folding of the backbone, optimizing values of  $ \phi $and $ \psi $and $ \omega $.
2. "Packing": Packing of the side chains of residues, optimizing values of $ \chi $etc.
3. "Designing": Designing of the protein sequence, optimizing the total energy by changing the identities of amino acids.   
(e.g. substitute an Arg with a Gly).

Monte Carlo simulation is the most basic and classic folding algorithm. 

The basic process of the Monte Carlo simulation is shared by all folding, packing and design process.

<details class="lake-collapse"><summary id="u85986045"><span class="ne-text">Basic Folding Process using Monte Carlo simulation</span></summary><p id="u1cf4ba01" class="ne-p"><a href="https://www.yuque.com/yuqueyonghu2r84jv/zmy8pd/gfkfxof6pygtuwkr" data-href="https://www.yuque.com/yuqueyonghu2r84jv/zmy8pd/gfkfxof6pygtuwkr" class="ne-link"><span class="ne-text">Note | PyRosetta Basics (6): Mover</span></a></p><h3 id="e64YK"><span class="ne-text" style="color: rgba(0, 0, 0, 0.87)">Making a trial move</span></h3><p id="u5e5dbfd4" class="ne-p"><span class="ne-text" style="font-size: 14px">Apply a </span><em><span class="ne-text" style="font-size: 14px">move</span></em><span class="ne-text" style="font-size: 14px"> (change) to the </span><em><span class="ne-text" style="font-size: 14px">pose</span></em><span class="ne-text" style="font-size: 14px"> (structure).</span></p><h3 id="irwtJ"><span class="ne-text">Scoring the move</span></h3><p id="u44f841d4" class="ne-p"><span class="ne-text" style="font-size: 14px">Compare the energy scores before and after the move. </span></p><p id="u9e89c1ac" class="ne-p"><span class="ne-text" style="font-size: 14px">(Therefore, the pose should be cloned as a backup before the move.)</span></p><p id="u24599695" class="ne-p"><span class="ne-text" style="font-size: 14px"></span></p><p id="u8cc42be0" class="ne-p"><span class="ne-text" style="font-size: 14px">And, compare the energy score after the </span><em><span class="ne-text" style="font-size: 14px">move</span></em><span class="ne-text" style="font-size: 14px"> and the lowest energy score in history. </span></p><p id="u020e8d13" class="ne-p"><span class="ne-text" style="font-size: 14px">If the score after this </span><em><span class="ne-text" style="font-size: 14px">move</span></em><span class="ne-text" style="font-size: 14px"> is the lowest, store the </span><em><span class="ne-text" style="font-size: 14px">pose </span></em><span class="ne-text" style="font-size: 14px">after this</span><em><span class="ne-text" style="font-size: 14px"> move </span></em><span class="ne-text" style="font-size: 14px">as the best pose, and renew the lowest energy score in history.</span></p><h3 id="TAn8b"><span class="ne-text">Desiding whether or not to accept the move: </span><strong><span class="ne-text">Metropolis criterion</span></strong></h3><p id="ucb6377ba" class="ne-p"><span class="ne-text">For the</span><strong><span class="ne-text"> </span></strong><span class="ne-text">decision step, we need to make a subroutine that either accepts or rejects the new conformatuon based on the </span><strong><span class="ne-text">Metropolis criterion</span></strong><span class="ne-text">. </span></p><p id="u555ca55a" class="ne-p"><span class="ne-text">The difference between the energy scores after (t+1) and before (t) a move is </span><span id="YTTPK" class="ne-math" style="font-size: 14px"><img src="https://cdn.nlark.com/yuque/__latex/140cea98d95645c2586ed9f29e31675b.svg"></span><span class="ne-text">.</span></p><p id="u4257265d" class="ne-p"><span class="ne-text">The Metropolis criterion has a probability of accepting a move as </span><span id="GuL5z" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/862838486b7050697b0dee3c6ecb5d6d.svg"></span><span class="ne-text">.</span></p><p id="ub1003231" class="ne-p"><span class="ne-text">When </span><span id="xfX93" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/fb3b0300200beee945a7c33ee9da1833.svg"></span><span class="ne-text">, the Metropolis criterion probability of accepting the move is </span><span id="BH3hC" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/862838486b7050697b0dee3c6ecb5d6d.svg"></span><span class="ne-text">.</span></p><p id="u03854502" class="ne-p"><span class="ne-text">When </span><span id="jAP7B" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/16d5d059d93f960fa208a43278bc6094.svg"></span><span class="ne-text">, the Metropolis criterion probability of accepting the move is </span><span id="U1a3A" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/54928c26406c51ef50311d6f0d99f3a6.svg"></span><span class="ne-text">. </span></p><p id="ucd6574c7" class="ne-p"><span class="ne-text">Use </span><span id="lfd2T" class="ne-math"><img src="https://cdn.nlark.com/yuque/__latex/825509854c41e934535f33223f7981f1.svg"></span><span class="ne-text"> . </span></p><p id="u05b7dc94" class="ne-p"><span class="ne-text"></span></p><p id="uf6ba059c" class="ne-p"><span class="ne-text">If the move is accepted, the </span><em><span class="ne-text">pose</span></em><span class="ne-text"> after the </span><em><span class="ne-text">move</span></em><span class="ne-text"> will be used for the next round of trial </span><em><span class="ne-text">move.</span></em></p><p id="ud6027a74" class="ne-p"><span class="ne-text">If not, the </span><em><span class="ne-text">pose</span></em><span class="ne-text"> before the </span><em><span class="ne-text">move </span></em><span class="ne-text">will be used for the next round of trial </span><em><span class="ne-text">move</span></em><span class="ne-text">.</span></p><h3 id="SJZe5"><span class="ne-text">Iterations</span></h3><p id="u892da1c2" class="ne-p"><span class="ne-text">For each iteration, the 3 steps above are excuted. </span></p><p id="u6a79e2e1" class="ne-p"><span class="ne-text">The final output of this program should be the </span><strong><span class="ne-text">lowest energy conformation</span></strong><span class="ne-text"> that is achieved </span><strong><span class="ne-text">at any point</span></strong><span class="ne-text"> during the simulation. </span></p><p id="u53388aed" class="ne-p"><span class="ne-text" style="color: #DF2A3F">NOTICE</span><span class="ne-text">: Search &quot;greedy algorithm&quot; and &quot;Monte Carlo&quot; for more information about the algorithm design. </span></p></details>
### Simulated annealing
Simulated anealing refer to the ramping of temperature `kT`in Python looping.

1. Set a high initial temperature, so that more barriers can be easily overcome.
2. By gradually decreasing the temperarture after each Python loop, the barriers becomes harder to overcome, and the search range gradually narrows down to a local minimum.
3. In this process, after a period of "heated" (active) searching across the energy landscape, the energy of a pose is expected to be gradually trapped in a deep trough (valley), which is a reliable local energy minimum (a possible global minimum).

### Simulated relaxing
Simulated relaxing refer to the ramping of energy score `E`(and the differences between scores `delta_E`) in Python looping.

1. Set a small initial weight of an energy term, so that the free energy (and the barriers) are underestimated, and thus the energy barriers are easier to be overcome.
2. By gradually increasing the weight after each Python loop, the energy estimation becomes more realistic, and the search range gradually narrows down to a local minimum.
3. In this process, after a period of the "relaxed" (loose) evaluation, energy of a pose is expected to be gradually trapped in a deep trough (valley), which is a reliable local energy minimum (a possible global minimum).

## Introduction to Packing
"Packing" is the process of optimizing the conformation of side chains (values of $ \chi_{n} $and other angles).

The "_packer_" is a subset of _movers _that apply packing operations to a pose_. _

We use "task operations" to configure the behaviors of _packers_:

+ The "task operations" are provided to a `TaskFactory`object.
+ The  `TaskFactory` is then passed to a _packer_.

## `init()`：Initialize From Rosetta Commandline
We can change some defalt Rosetta (C++) commandline settings through

```python
pyrosetta.rosetta.init()
```

E.g.

```python
# we can change some defalt settings through
# pyrosetta.rosetta.init()
# or, write as:
# pyrosetta.init()
# For example
init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false -relax:default_repeats 5 -no_fconfig')

'''
In this initialization, we
1. input the scheme for antibody nomenclature;
2. use the PDB ligand definitions to load all standard residue descriptions;
3. keep the relax mode as `FastRelax` (which is default);
4. and set the `FastRelax` repeat times to 5 (by default it is already 5, can be customized);
5. and do not load the [common] config file if present.
'''
```

## `TaskFacrtory`class: prototype of `PackerTask`
### Introduction
Instances of`PackerTask`class are used to determine the specific tasks in packing.

`TaskFactory`is the prototype of `PackerTask` class.

Every time the protein is packed for a single round, the `TaskFactory` will generate what is called the `PackerTask`.



We do NOT use `PackerTask` directly.

We can set `TaskOperations`to the `TaskFactory`, and all `PackerTask`generated by this`TaskFactory` will be controlled by these `TaskOpterations`.

So bascically a `TaskFactory`object can be seen as "a list of task operations".

Some `TaskOperations`can respond to changes in the pose, so we do not set `TaskOperations` directly to `PackerTask`in each packing step.

Instead, we use `TaskFactory`to dynamically generate `PackerTask`objects in each round, behind the scene. 

### Instantiation & Configuration(`push_back`)
```python
my_tf = TaskFactory() # Create an empty task factory

# And set task operations to the task factory:
my_tf.push_back(operation.InitializeFromCommandline()) 
# IMPORTANT! Use configurations declared in `init()`
my_tf.push_back(operation.RestrictToRepacking()) 
# Disable "designing of side chains" so that the sequence will be untouched
```

### Regional Packing: `selections` & `PreventRepackingRLT`
We can also use `selections`for regional packing

E.g. Selecting CDR H1 region of an antibody

```python
# Selecting CDR H1 region of an antibody
cdr_selector = CDRResidueSelector() # Instantiation
cdr_selector.set_cdr(h1) # Select the H1 loop only
# the namespace of `h1` was imported by commandlines that were passed to init() 

# Selecting the Neiborhood (in 3D space) residue of CDR H1
nbr_selector = selections.NeighborhoodResidueSelector() # Instantiation
nbr_selector.set_focus_selector(cdr_selector) # Find neighbors of cdr_selector (which is the H1 loop)
nbr_selector.set_include_focus_in_subset(True) # Include the H1 loop in `nbr_selector` as well

#"RLT": Residue Level Task Operation
prevent_repacking_rlt = operation.PreventRepackingRLT() # Instantication

prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True)
#`True` indicates here that we are flipping the selection.  
# So that we are turning off everything but the CDR and its neighbors.

# Assume that we have a `TaskFactory` named `my_tf`
my_tf.push_back(prevent_subset_repacking)
```

### `TaskFacrtory` Methods
#### Add operations: `TaskFacrtory.push_back`
```python
my_tf.push_back(my_operation) 
# E.g.
my_tf.push_back(operation.InitializeFromCommandline()) 
```

#### Clear all operations: `TaskFacrtory.clear`
```python
my_tf.clear()
```

### Resources
[https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.pack.task.operation.html](https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.pack.task.operation.html)



## `PackRotamersMover`class
A _packer_ is a _mover_ that carries out packing under the instruction of `PackerTask`.

### Instantiation & Configurations
```python
# Assume that we have a `TaskFactory` object called `my_tf`
my_packer = pack_min.PackRotamersMover() # Create a `PackRotamersMover` object
my_packer.task_factory(my_tf) # Pass the `TaskFactory` object to the packer

#Note that we are not passing a scorefunction here.  We will use the default, cmd-line scorefunction, 
# which can be accessed through `rosetta.core.scoring.get_score_function()`
# We use a scorefunction later. 
```

### Apply the packer to a pose
```python
#Run the packer. (Note this may take a few minutes)
my_packer.apply(pose)
```

