#!/usr/bin/env python3

import sys
from ase.io import read, write
from ase.visualize import view
from acat.adsorption_sites import ClusterAdsorptionSites
from acat.adsorbate_coverage import ClusterAdsorbateCoverage
from acat.build import add_adsorbate_to_site, remove_adsorbates_from_sites

filename = '../CONTCAR'
atoms = read(filename)
print(atoms)

atoms.pbc = False
cas = ClusterAdsorptionSites(atoms, surrogate_metal='Cu', tol=0.1) #, composition_effect=True) #, label_sites=True)

# unique sites
unique_sites = cas.get_unique_sites() #unique_composition=True)
print(unique_sites)

for sidx, site in enumerate(unique_sites):
  add_adsorbate_to_site(atoms, adsorbate='N', site=site)
  #print(f"Adsorbing on site: {sidx+1:03}", site['site'], site['composition'], site['label'])
  print(f"Adsorbing on site: {sidx+1:03}", site['site'], site['label'])
  #filename_ads = f"adsorption_{sidx+1:03}"+"-"+site['site']+"-"+site['composition']+".xyz"
  filename_poscar = f"POSCAR_{sidx+1:03}"
  atoms.write(filename_poscar, format='vasp', vasp5=True)  
  cac = ClusterAdsorbateCoverage(atoms, cas)
  occupied_sites = cac.get_sites() #occupied_only=True)
  remove_adsorbates_from_sites(atoms, sites=occupied_sites)
