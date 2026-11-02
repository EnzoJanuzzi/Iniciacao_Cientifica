#!/usr/bin/env python3
from ase.io import read, write
from ase.constraints import FixAtoms

# 1. Read your optimized structure
atoms = read('POSCAR')

# 2. Identify which atoms to freeze (Cu and Ru)
frozen_indices = [atom.index for atom in atoms if atom.symbol in ['Cu', 'Ru']]

# 3. Apply the constraint (F F F)
constraint = FixAtoms(indices=frozen_indices)
atoms.set_constraint(constraint)

# 4. Save the new POSCAR with Selective Dynamics enabled
write('POSCAR_freq', atoms, format='vasp', vasp5=True)

print(f"Success! Froze {len(frozen_indices)} metal atoms.")
print(f"Left {len(atoms) - len(frozen_indices)} adsorbate atoms free to vibrate.")
