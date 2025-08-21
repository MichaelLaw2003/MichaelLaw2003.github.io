title: 'PyRosetta Basics (6): Movers for Backbone Folding'
tags:
  - PyRosetta
  - biochemistry
  - protein design
  - note
  - protein
categories:
  - protein design
---
# PyRosetta Basics (6): Movers for Backbone Folding
```python
!pip install pyrosettacolabsetup
import pyrosettacolabsetup; pyrosettacolabsetup.install_pyrosetta()
import pyrosetta; pyrosetta.init()
from pyrosetta import *
from pyrosetta.teaching import *
init()
```

## Introduction
One of the most basic operations in protein structure and design algorithms is manipulation of the protein conformation. 

In Rosetta, these manipulations are organized into `mover`s. 

A `Mover` object simply changes the conformation of a given pose. 

It can be simple, like a single φ or ψ angle change, or complex, like an entire refinement protocol.

[Movers (RosettaScripts)](https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/Movers-RosettaScripts)

## `MoveMap`class
<font style="color:rgba(0, 0, 0, 0.87);">Most </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">Movers</font>`<font style="color:rgba(0, 0, 0, 0.87);"> require a </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">MoveMap</font>`<font style="color:rgba(0, 0, 0, 0.87);"> object to specify which degrees of freedom are fixed and which are free to change. </font>

<font style="color:rgba(0, 0, 0, 0.87);">Example1:</font>

```python
my_movemap = MoveMap()
my_movemap.set_bb(True) # Allowing a mover to change the backbone structure
```

Example2:

```python
movemap.set_bb(False) # Prohibit a mover to change the backbone structure
movemap.set_bb(50, True) # Allow a mover to chage the backbone structure of residue #50
movemap.set_bb(51, True) # Allow a mover to chage the backbone structure of residue #51
# Only the backbone structure at residue #50 and #51 can be modified by a mover
```

## Basic Backbone Folding Process
### <font style="color:rgba(0, 0, 0, 0.87);">Making a trial move</font>
Apply a _move_ (change) to the _pose_ (structure).

### Scoring the move
Compare the energy scores before and after the move. 

(Therefore, the pose should be cloned as a backup before the move.)



And, compare the energy score after the _move_ and the lowest energy score in history. 

If the score after this _move_ is the lowest, store the _pose _after this_ move _as the best pose, and renew the lowest energy score in history.

### Desiding whether or not to accept the move
For the** **decision step, we need to make a subroutine that either accepts or rejects the new conformatuon based on the **Metropolis criterion**. 



The difference between the energy scores after (t+1) and before (t) a move is $ \Delta E = E_{t+1} - E_{t} $.

The Metropolis criterion has a probability of accepting a move as $ P = \exp( -\Delta E / kT ) $.

When $ ΔE ≥ 0 $, the Metropolis criterion probability of accepting the move is $ P = \exp( -\Delta E / kT ) $.

When $ ΔE < 0 $, the Metropolis criterion probability of accepting the move is $ P = 1 $. 

Use $ kT = 1 \text{\ Rosetta \ Energy \ Unit \ (REU)} $ . 



<font style="color:#DF2A3F;">NOTICE</font>: Search "Boltzmann function" and "Simulated annealing" for more information about $ kT $.

If the move is accepted, the _pose_ after the _move_ will be used for the next round of trial _move._

If not, the _pose_ before the _move _will be used for the next round of trial _move_.

### Iterations
For each iteration, the 3 steps above are excuted. 

The final output of this program should be the **lowest energy conformation** that is achieved **at any point** during the simulation. 

<font style="color:#DF2A3F;">NOTICE</font>: Search "greedy algorithm" and "Monte Carlo" for more information about the algorithm design. 

## Basic backbone movers
### `SmallMover` & `ShearMover` classes
| **Mover** | **Description** |
| --- | --- |
| `SmallMover` | Makes "small-move-style" torsion moves (no propagation minimization) |
| `ShearMover` | Makes "shear-style" torsion moves that **minimize downstream propagation** |


<font style="color:rgba(0, 0, 0, 0.87);">For convenience, the </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">SmallMover</font>`<font style="color:rgba(0, 0, 0, 0.87);"> and </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">ShearMover</font>`<font style="color:rgba(0, 0, 0, 0.87);"> can do multiple rounds of perturbation. </font>

<font style="color:rgba(0, 0, 0, 0.87);">They also check that the new φ/ψ combinations are within an allowable region of the </font>**<font style="color:rgba(0, 0, 0, 0.87);">Ramachandran plot</font>**<font style="color:rgba(0, 0, 0, 0.87);"> by using a </font>**<font style="color:rgba(0, 0, 0, 0.87);">Metropolis acceptance criterion</font>**<font style="color:rgba(0, 0, 0, 0.87);"> based on the rama score component change. (The </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">rama</font>`<font style="color:rgba(0, 0, 0, 0.87);"> score is a statistical score from Simons et al. 1999, parametrized by bins of φ/ψ space.) </font>

<font style="color:rgba(0, 0, 0, 0.87);">Because they use the Metropolis criterion, we must also supply </font>`<font style="color:rgba(0, 0, 0, 0.87);">kT</font>`<font style="color:rgba(0, 0, 0, 0.87);">.</font>

```python
kT = 1.0 # temperature
n_moves = 1 # iterations

movemap = MoveMap()
movemap.set_bb(True)

small_mover = SmallMover(movemap, kT, n_moves)
shear_mover = ShearMover(movemap, kT, n_moves)
```

### `BackrubMover`
| **<font style="color:rgb(51, 51, 51);">Mover</font>** | **<font style="color:rgb(51, 51, 51);">Description</font>** |
| --- | --- |
| <font style="color:#000000;">BackrubMover</font> | <font style="color:rgb(51, 51, 51);">Makes local rotations around two backbone atoms</font> |


### <font style="color:rgb(51, 51, 51);">Minimization Mover </font>`<font style="color:rgb(51, 51, 51);">MinMover</font>`class
The `MinMover` carries out a **gradient-based** minimization to find the **nearest local minimum** in the energy function, such as that used in one step of the Monte-Carlo-plus-Minimization algorithm of Li & Scheraga. 

```python
min_mover = MinMover()
```

## `MonteCarlo`class
### Introduction
<font style="color:rgba(0, 0, 0, 0.87);">The </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">MonteCarlo</font>`<font style="color:rgba(0, 0, 0, 0.87);"> object is an encapsulated object that creates a whole MonteCarlo simulation.</font>

<font style="color:rgba(0, 0, 0, 0.87);">That is, it can decide whether to accept or reject a trial conformation, and it keeps track of the lowest-energy conformation and other statistics about the search. </font>

<font style="color:rgba(0, 0, 0, 0.87);">Having the Monte Carlo operations packaged together is convenient, especially if we want multiple Monte Carlo loops to nest within each other or to operate on different parts of the protein.</font>

```python
# Assume that we have instantiated a mover `my_mover`

# Instantiate a `MonteCarlo` object `mc`
mc = MonteCarlo(my_pose, scorefxn, kT)
# my_pose: a Pose object
# scorefxn: a ScoreFunction object
# kT: a parameter - temperature; int or float 

# apply the mover
my_mover.apply(my_pose) 

# deside whether or not to accept the move
flag = mc.boltzmann(my_pose)
# flag is boolean
# True: accept
# False: reject

# In practice, these 3 steps should be written into a loop.
```



A Monte Carlo Simulation process is (in pseudo code):

```python
def monte_carlo_sampling(initial_state, kT, N_iteration):
    current_state = initial state
    for i in range(N_iteration):
        new_state = random_move(current_state)
        #fe_delta: delta (difference of) free energies
        fe_delta = free_energy(new_state)-free_energy(current_state) 
        if fe_delta < 0 or uniform(1) < exp(-fe_delta / kT): 
            # exp(-fe_delta / kT): boltzmann distribution
            current_state = new_state
```

### `MonteCarlo` Methods
```python
mc.boltzmann(my_pose) # See above
mc.show_scores()
mc.show_counters()
mc.show_state()
```

### `TrialMover`class
<font style="color:rgba(0, 0, 0, 0.87);">A </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">TrialMover</font>`<font style="color:rgba(0, 0, 0, 0.87);"> combines a specified </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">Mover</font>`<font style="color:rgba(0, 0, 0, 0.87);"> with a </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">MonteCarlo</font>`<font style="color:rgba(0, 0, 0, 0.87);"> object.</font>

<font style="color:rgba(0, 0, 0, 0.87);">Each time a </font>`<font style="color:rgba(0, 0, 0, 0.87);background-color:rgb(238, 238, 238);">TrialMover</font>`<font style="color:rgba(0, 0, 0, 0.87);"> is called, it performs a trial move and tests that move' s acceptance with the </font>`<font style="color:rgba(0, 0, 0, 0.87);">MonteCarlo</font>`<font style="color:rgba(0, 0, 0, 0.87);"> object. </font>

<font style="color:rgba(0, 0, 0, 0.87);">It is designed to test the effects of a </font>`<font style="color:rgba(0, 0, 0, 0.87);">Mover</font>`<font style="color:rgba(0, 0, 0, 0.87);">.</font>

```python
trial_mover = TrialMover(small_mover, mc)

for i in range(10):
    trial_mover.apply(test)
    
print(trial_mover.num_accepts())
print(trial_mover.acceptance_rate())

# After the trial, information about the trial can also be visited by
# mc.show_state()
# so that different movers can be compared
```



## `<font style="color:rgba(0, 0, 0, 0.87);">SequenceMover</font>`class and `RepeatMover`class
```python
seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)
```

```python
n_repeats = 3
repeat_mover = RepeatMover(my_mover, n_repeats)
```

