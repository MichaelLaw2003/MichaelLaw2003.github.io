title: 'PyRosetta Basics (1): Pose'
date: '2025-08-21 12:42:40'
updated: '2025-08-21 12:42:42'
tags:
  - PyRosetta
  - protein design
  - biochemistry
  - protein
  - note
categories:
  - protein design
---
# PyRosetta Basics (1): Pose
```python
from pyrosetta import *
init()
#import pyrosetta
#pyrosetta.init()
```

## `Pose`class
### Introduction
1. Simplified explanation: 
    1. `Pose`(`pyrosetta.Pose`) is defined as a Python class
    2. `Pose`is the core concept of PyRosetta
    3. An instance of the`Pose`class is refered to as a "pose object".
    4. A pose object is a protein structure.
2. Further explanation:
    1. A pose object may also include informations about related ions, ligand and other molecules.
    2. A pose object may have more than one protein chains (quaternary protein, or protein complex)
    3. Like any other Python classes, the`Pose`class has methods and attributes.
    4. Class methods can be visited from the instances of the class.

### `Pose` instantiation
#### Initiate an instance of `Pose` from a PDB file
```python
pose_1 = pose_from_pdb("./path_to_pdb/5tj3.pdb") 
# pyrosetta.pose_from_pdb
# load '5tj3.pdb' from file
```

#### Initiate an instance of `Pose` from [https://www.rcsb.org/](https://www.rcsb.org/)
```python
pose_2 = pose_from_rcsb("5TJ3")
# pyrosetta.pose_from_rcsb
# load 5TJ3 from online database
```

#### Initiate an instance of `Pose`from a string (a sequence)
```python
pose_3 = pose_from_sequence("AAAAAAAAAA")
# pyrosetta.pose_from_sequence
# load a "poly alanine sequence" from string
```

#### Initiate an empty instance of `Pose` (And create a clone of a `pose` object)
```python
pose_4 = Pose() # Now pose_4 is empty

# we can now use pose_4 to create a deep clone (copy) of pose_1
pose_4.assign(pose_1)
```

### `Pose`methods
#### Copy a `pose`object via `Pose.assign`
```python
# Assume that we have a Pose instance called `my_pose`
my_pose_clone = Pose() # Now my_pose_clone is empty

#create a deep clone (copy) of my_pose
my_pose_clone.assign(my_pose)


```

<font style="color:#DF2A3F;">WARNING</font>: Avoid writing something like:

```python
my_pose_clone = my_pose 
# This is incorrect
```

Because in this way, my_pose_clone is pointed to value of my_pose, and modification to my_pose_clone will change my_pose.

#### Copy a `pose`object via `Pose.clone`
```python
# assume that we have a `pose` object called `my_pose`
my_pose_clone = my_pose.clone()
```

#### Get sequence
```python
# `my_pose` is an instance of Pose
my_seq = my_pose.sequence()
# `my_seq` is a string, storing the sequence of my_pose
```

#### Get annotated sequence
```python
my_seq_annotated = my_pose.annotated_sequence()
# `my_seq` is a string, storing the annotated sequence of my_pose
# annotations are like: '...T[THR:phosphorylated]...'
```

#### Get total residue count
```python
my_total_residue_count = my_pose.total_residue()
# `my_total_residue_count` is an integer
```

#### Get dihedral angles
```python
# `residue_id` is the ID (in pose) of a residue, type is integer
my_pose.phi(residue_id) # C-CaN-C dihedral
# >>> −64.8
my_pose.psi(residue_id) # N-CaC-N dihedral
# >>> -41.0
my_pose.chi(1, residue_id) # 1st C-C dihedral in side chain
# >>> -82.8
```

#### Get an `Energies`object
```python
my_energies = my_pose.energies()

print(my_energies.show(24)) # Print all energies terms (unweighted) of residue 24
```

#### Get a`Residue` object
```python
residue_20th = my_pose.residue(20)
# `residue_20th` is an object that stores the 20th residue of `my_pose`
```

See "Residue objects" in this note 

<font style="color:#DF2A3F;">WARNING</font>: "the 20th residue in pose" is not equal to "the 20th residue in PDB", because a pdb file may contain more than one chain.

#### Get a`PdbInfo`object
This object is a bridge between a `Pose` object and a PDB file:

1. This object stores information from the original PDB file (if the pose is from a pdb file or from the pdb database).
2. Or, this object stores information that can be written into a PDB file.

```python
my_pdb_info = my_pose.pdb_info()
```

See "PdbInfo objects" in this note 

#### Get a`Conformation`object
```python
my_conformation = my_pose.conformation() # Get conformation object
```

#### Set dihedral angles (a basic mover)
```python
my_pose.set_phi(24, -64.8) # set phi of residue #24 to -64.8
my_pose.set_psi(24, -41.0) # set psi of residue #24 to -41.0
```

#### Set foldtree (requires a pre-defined `FoldTree`object)
```python
# Assume that we have a `FoldTree` object called `my_ft`
my_pose.fold_tree(my_ft)
```

#### Output (Save) a pose as a PDB file
```python
my_pose.dump_pdb('/outputs/my_pose_arxiv.pdb')
```

## `Residue` class
`Residue`is a secondary class of `Pose`.

A `residue`object stores information about a residue.

<font style="color:#DF2A3F;">NOTICE</font>: An ion, or a small molecule ligand, etc. is also seen as a residue.

### Get residue name
```python
residue_20th = my_pose.residue(20)
# my_pose is a `Pose` object
# residue_20th is a `Residue` object
residue_20th_name = residue_20th.name()
# residue_20th_name is a string, storing the uppercased 3-letter name
# e.g. 'ASP'
```

### Get residue property: `is_xxx`booleans
```python
residue_24th = my_pose.residue(24) # Get #24 residue from my_pose
residue_24th_is_charged = residue_24th.is_charged() # Boolean, `True` if it is charged
residue_24th_is_aa = residue_24th.is_protein() # Boolean, `True` if it is an amino acid
```

### Get atom index
```python
residue_24th = my_pose.residue(24) # get #24 residue
carbon_alpha_of_residue_24th = residue_24th.atom_index("CA")
```

`"N"`is the nitrogen of main-chain amino group.

`"CA"`is the central carbon.

`"C"`is the carbon of main-chain carboxyl group.

For atom nomenclatures, search for "amino acid structure".

### Get atom xyz coordinates (for vector-related calculations)
```python
# `residue_24th` is a residue in `my_pose`
N_xyz = residue_24th.xyz("N")
# N_xyz is a vector (length = 3)
```



## `PdbInfo`class
`PdbInfo`is a secondary class of `Pose`.

This class is a bridge between a `Pose` object and a PDB file:

1. A `pdb_info`object stores information from the original PDB file (if the pose is from a pdb file or from the pdb database).
2. Or, this object stores information that can be written into a PDB file.

### Convert residue ID (number): `pdb2pose` and `pose2pdb`
The ID (number) of a residue starts at #1 in pose.

The ID of a residue in a `pose`object is different from its ID in a `.pdb`file.

```python
# In PDB file, a residue of interect is in "chain A", and its ID (number) in the PDB file is 24
# But its ID in the `Pose` object (my_pose) is unknown.
my_pdb_info = my_pose.pdb_info()
my_residue_id_in_pose = my_pdb_info.pdb2pose('A', 24) # get the residue ID in pose, stored as an integer
my_residue = my_pose.residue(my_residue_id_in_pose) # get the `Residue` object by its ID in pose

# alternatively, we can call the methods in a chain: 
# my_residue_id_in_pose = my_pose.pdb_info().pdb2pose('A', 24)
```

```python
# On the other hand, we can get the corresponding PDB chain infomation and ID information for a specific residue
# In a `Pose` object, a residue is at #1 position (its ID is 1), but its location in corresponding PDB file is unknown
my_residue_info_in_pdb = my_pdb_info.pose2pdb(1)
# my_residue_info_in_pdb is a string that looks like '24 A'
```

### Get chain information & Get number (ID) information
```python
# a `Pose` instance called `my_pose`, and a residue of interest whose ID (in pose) is 1
my_pose.pdb_info().chain(1)
# >>> A
# The #1 residue is in 'chain A'
my_pose.pdb_info().number(1)
# >>> 24
# The ID of #1 (in pose) residue in corresponding PDB file is 24
```

### Set `PdbInfo`
```python
# assume that my_pose is a `pose` object
my_pose.pdb_info().name("test") 
# set the name of this pose in its pdb_info as "test"
# if it is sent to PyMOL, it will be displayed as "test" in the object list
```

## `Conformation`class
### Get bond length & Get bond angle
```python
# given a `residue_id` (ID of a residue):
residue_28th = my_pose.residue(residue_id) # Get residue 
N28 = AtomID(residue_28th.atom_index("N"), residue_id) # Construct an `AtomID` object for backbone nitrogen
CA28 = AtomID(residue_28th.atom_index("CA"), residue_id)
C28 = AtomID(res_28.atom_index("C"), resid)

my_conformation = my_pose.conformation() # Get conformation object

N_CA_28_bond_length = my_conformation.bond_length(N28, CA28) # get bond length; type is float
bb_angle = my_conformation.bond_angle(N28, CA28, C28) # get bond angle; type is float; in degree
```

