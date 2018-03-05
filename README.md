# ProtonHopping
Data and code for Proton Hopping project. The gaol of this project is to develop a method to calculate
energies of proton hopping in hydrogen-bonded water-water and protein-water systems. See a brief project summary [here](http://kamranhaider.org/projects.html)

# Current Status
MCCE simulation protocol under optimization. Current inputs can be reproduced as follows. All scripts should be run from
`scripts` subfolder.

### Step 1:
Generate MCCE compatible structure of gramicidin channel.
```python
python prep_gramicidin_structure.py
``` 
### Step 2:
Insert gramicidin into a membrane (using IPECE program).
```python
python insert_membrane.py
``` 
### Step 3:
Re-build full gramicidn-membrane MCCE compatible structure
```python
python prep_mem_gramicidin_structure.py
``` 
### Step 4:
Fix membrane atom names so that they conform to MCCE atom typing for membrane atoms.
```python
python fix_mem_atom_names.py
``` 
### Step 5:
Generate a rectangular gird and place water oxygen atoms, remove oxygens that
are in vdw clash with the system.
```python
python place_water_oxygens.py
```
### Step 6:
Generate water geometries from the oxygen atoms using tleap.
```bash
tleap -f tleap_type_waters.in
```
### Step 7:
For each oxygen atom, enumerate orientational conformers, resulting in N conformer per water,
where N can be 5, 10, 15, 20, 25.
```python
python generate_water_conformers.py
``` 



