#!/usr/bin/env python3
"""
Comprehensive Feature Extraction Script for Adsorbate-Nanoparticle Calculations.

Compatible with VASP outputs (CONTCAR, OUTCAR, vasprun.xml, OSZICAR).
Extracts:
  - Identification: intermediate, pathway, stage (screening/screening_final), conformation
  - Energetics: total energy (eV), energy per atom, energy per metal atom, Fermi energy (eV)
  - Convergence: electronic/ionic convergence, number of ionic steps, max force (eV/A)
  - Structure & Geometry: NP radius of gyration, NP center of mass, surface/core radii
  - Adsorbate Features: formula, atomic coordinates (x,y,z), center of mass, internal bond lengths,
    binding atom, nearest metal element, binding distance, coordination number, adsorption site (top/bridge/hollow)
  - Magnetization: total, absolute, core (Cu) mean, surface (Ru) mean, adsorbate total, binding metal atom
  - Electronic DOS: d-band center and d-band width for Ru (surface), Cu (core), and total metal.
"""

import os
import sys
import glob
import json
import argparse
import numpy as np
import pandas as pd
from ase.io import read
from pymatgen.io.vasp import Outcar, Vasprun
from pymatgen.core import Element
from pymatgen.electronic_structure.core import OrbitalType, Spin

METAL_ELEMENTS = {'Cu', 'Ru', 'Pt', 'Pd', 'Au', 'Ag', 'Fe', 'Co', 'Ni', 'Rh', 'Ir'}

def compute_d_band_stats(cdos, element_symbol, efermi):
    """
    Computes the d-band center and d-band width for a given element relative to Fermi level.
    Integrates occupied states (E <= efermi, i.e., E - efermi <= 0).
    """
    try:
        spd = cdos.get_element_spd_dos(Element(element_symbol))
        if OrbitalType.d not in spd:
            return np.nan, np.nan
        d_dos = spd[OrbitalType.d]
        energies = d_dos.energies - efermi

        densities = np.array(d_dos.densities[Spin.up])
        if Spin.down in d_dos.densities:
            densities = densities + np.array(d_dos.densities[Spin.down])

        mask = energies <= 0
        e_occ = energies[mask]
        rho_occ = densities[mask]

        norm = np.trapz(rho_occ, e_occ)
        if norm <= 1e-8:
            return np.nan, np.nan

        d_center = np.trapz(e_occ * rho_occ, e_occ) / norm
        d_width = np.sqrt(np.trapz((e_occ - d_center)**2 * rho_occ, e_occ) / norm)
        return float(d_center), float(d_width)
    except Exception:
        return np.nan, np.nan


def extract_features_from_calc(calc_dir):
    """
    Extracts features from a single calculation directory.
    """
    contcar_path = os.path.join(calc_dir, "CONTCAR")
    outcar_path = os.path.join(calc_dir, "OUTCAR")
    vasprun_path = os.path.join(calc_dir, "vasprun.xml")

    if not os.path.exists(contcar_path) or os.path.getsize(contcar_path) == 0:
        contcar_path = os.path.join(calc_dir, "POSCAR")
        if not os.path.exists(contcar_path):
            return None

    if not os.path.exists(outcar_path):
        return None

    abs_path = os.path.abspath(calc_dir)
    path_parts = abs_path.split(os.sep)
    conf_name = os.path.basename(abs_path)

    stage = "unknown"
    intermediate = "unknown"
    pathway = "unknown"
    try:
        conf_idx = path_parts.index(conf_name)
        stage = path_parts[conf_idx - 1]
        if path_parts[conf_idx - 2] == "dock":
            intermediate = path_parts[conf_idx - 3]
            pathway = path_parts[conf_idx - 4]
        else:
            intermediate = path_parts[conf_idx - 2]
    except Exception:
        pass

    try:
        atoms = read(contcar_path)
    except Exception:
        return None

    n_total_atoms = len(atoms)
    symbols = atoms.get_chemical_symbols()
    positions = atoms.positions

    metal_indices = [i for i, s in enumerate(symbols) if s in METAL_ELEMENTS]
    ads_indices = [i for i, s in enumerate(symbols) if s not in METAL_ELEMENTS]

    n_metal = len(metal_indices)
    n_adsorbate = len(ads_indices)

    if n_metal > 0:
        metal_pos = positions[metal_indices]
        np_com = np.mean(metal_pos, axis=0)
        metal_dists_from_com = np.linalg.norm(metal_pos - np_com, axis=1)
        r_gyration_np = float(np.sqrt(np.mean(metal_dists_from_com**2)))

        cu_indices = [i for i in metal_indices if symbols[i] == 'Cu']
        ru_indices = [i for i in metal_indices if symbols[i] == 'Ru']

        core_mask_local = metal_dists_from_com < 3.0
        surf_mask_local = ~core_mask_local

        np_surf_min_r = float(np.min(metal_dists_from_com[surf_mask_local])) if any(surf_mask_local) else np.nan
        np_surf_mean_r = float(np.mean(metal_dists_from_com[surf_mask_local])) if any(surf_mask_local) else np.nan
        np_surf_max_r = float(np.max(metal_dists_from_com[surf_mask_local])) if any(surf_mask_local) else np.nan
    else:
        np_com = np.zeros(3)
        r_gyration_np = np.nan
        np_surf_min_r = np_surf_mean_r = np_surf_max_r = np.nan
        cu_indices = []
        ru_indices = []

    ads_formula = "".join([symbols[i] for i in ads_indices])
    ads_coords_dict = {}
    closest_metal_element = "None"
    closest_metal_idx = -1
    binding_distance = np.nan
    binding_ads_symbol = "None"
    coord_number_2_4A = 0
    ads_com = [np.nan, np.nan, np.nan]
    internal_bonds = {}

    if n_adsorbate > 0:
        ads_pos = positions[ads_indices]
        ads_com = np.mean(ads_pos, axis=0).tolist()

        for idx in ads_indices:
            sym = symbols[idx]
            pos = positions[idx]
            ads_coords_dict[f"{sym}_{idx}"] = [round(float(c), 4) for c in pos]

        for i, idx1 in enumerate(ads_indices):
            for j, idx2 in enumerate(ads_indices):
                if j > i:
                    d_int = float(np.linalg.norm(positions[idx1] - positions[idx2]))
                    if d_int < 2.2:
                        pair_name = f"bond_{symbols[idx1]}_{symbols[idx2]}"
                        internal_bonds[pair_name] = round(d_int, 4)

        min_overall_dist = 999.0
        for idx in ads_indices:
            d_metals = np.linalg.norm(positions[metal_indices] - positions[idx], axis=1)
            min_m_idx = np.argmin(d_metals)
            min_d = d_metals[min_m_idx]
            if min_d < min_overall_dist:
                min_overall_dist = min_d
                binding_ads_symbol = symbols[idx]
                closest_metal_idx = metal_indices[min_m_idx]
                closest_metal_element = symbols[closest_metal_idx]
                binding_distance = float(min_d)
                coord_number_2_4A = int(np.sum(d_metals <= 2.45))

    if coord_number_2_4A == 1:
        site_type = "top"
    elif coord_number_2_4A == 2:
        site_type = "bridge"
    elif coord_number_2_4A >= 3:
        site_type = "hollow"
    elif n_adsorbate == 0:
        site_type = "clean_np"
    else:
        site_type = "weak_or_physisorbed"

    outcar = Outcar(outcar_path)
    total_energy = float(outcar.final_energy) if outcar.final_energy is not None else np.nan
    energy_per_atom = total_energy / n_total_atoms if n_total_atoms > 0 else np.nan
    energy_per_metal_atom = total_energy / n_metal if n_metal > 0 else np.nan

    total_mag = float(outcar.total_mag) if outcar.total_mag is not None else 0.0
    site_mags = [m['tot'] for m in outcar.magnetization] if outcar.magnetization else [0.0] * n_total_atoms
    abs_mag = float(np.sum(np.abs(site_mags))) if site_mags else 0.0

    mean_cu_mag = float(np.mean([site_mags[i] for i in cu_indices])) if cu_indices else 0.0
    mean_ru_mag = float(np.mean([site_mags[i] for i in ru_indices])) if ru_indices else 0.0
    ads_tot_mag = float(np.sum([site_mags[i] for i in ads_indices])) if ads_indices else 0.0
    binding_metal_mag = float(site_mags[closest_metal_idx]) if closest_metal_idx >= 0 and closest_metal_idx < len(site_mags) else 0.0

    converged_electronic = None
    converged_ionic = None
    n_ionic_steps = 0
    max_force = np.nan
    efermi = np.nan
    d_center_Ru = np.nan
    d_width_Ru = np.nan
    d_center_Cu = np.nan
    d_width_Cu = np.nan

    if os.path.exists(vasprun_path) and os.path.getsize(vasprun_path) > 0:
        try:
            vasprun = Vasprun(vasprun_path, parse_dos=True)
            converged_electronic = vasprun.converged_electronic
            converged_ionic = vasprun.converged_ionic
            n_ionic_steps = len(vasprun.ionic_steps)
            if n_ionic_steps > 0:
                forces = vasprun.ionic_steps[-1].get('forces', [])
                if len(forces) > 0:
                    max_force = float(np.max(np.linalg.norm(forces, axis=1)))

            efermi = float(vasprun.efermi)
            cdos = vasprun.complete_dos
            if cdos is not None:
                d_center_Ru, d_width_Ru = compute_d_band_stats(cdos, 'Ru', efermi)
                d_center_Cu, d_width_Cu = compute_d_band_stats(cdos, 'Cu', efermi)
        except Exception:
            efermi = getattr(outcar, 'efermi', np.nan)

    record = {
        "pathway": pathway,
        "intermediate": intermediate,
        "stage": stage,
        "conformation": conf_name,
        "folder_path": abs_path,
        "total_energy_eV": total_energy,
        "energy_per_atom_eV": energy_per_atom,
        "energy_per_metal_eV": energy_per_metal_atom,
        "e_fermi_eV": efermi,
        "converged_electronic": converged_electronic,
        "converged_ionic": converged_ionic,
        "n_ionic_steps": n_ionic_steps,
        "max_force_eV_A": max_force,
        "n_total_atoms": n_total_atoms,
        "n_metal_atoms": n_metal,
        "r_gyration_metal_A": r_gyration_np,
        "np_com_x": round(float(np_com[0]), 4),
        "np_com_y": round(float(np_com[1]), 4),
        "np_com_z": round(float(np_com[2]), 4),
        "np_surf_min_r": np_surf_min_r,
        "np_surf_mean_r": np_surf_mean_r,
        "np_surf_max_r": np_surf_max_r,
        "adsorbate_formula": ads_formula,
        "n_adsorbate_atoms": n_adsorbate,
        "binding_adsorbate_atom": binding_ads_symbol,
        "closest_metal_element": closest_metal_element,
        "closest_metal_index": closest_metal_idx,
        "binding_distance_A": binding_distance,
        "coordination_number_2_4A": coord_number_2_4A,
        "adsorption_site": site_type,
        "ads_com_x": round(float(ads_com[0]), 4) if not np.isnan(ads_com[0]) else np.nan,
        "ads_com_y": round(float(ads_com[1]), 4) if not np.isnan(ads_com[1]) else np.nan,
        "ads_com_z": round(float(ads_com[2]), 4) if not np.isnan(ads_com[2]) else np.nan,
        "ads_internal_bonds": json.dumps(internal_bonds),
        "ads_coordinates": json.dumps(ads_coords_dict),
        "total_mag_muB": total_mag,
        "abs_mag_muB": abs_mag,
        "core_Cu_mag_mean": mean_cu_mag,
        "surf_Ru_mag_mean": mean_ru_mag,
        "binding_metal_mag": binding_metal_mag,
        "adsorbate_total_mag": ads_tot_mag,
        "d_band_center_Ru_eV": d_center_Ru,
        "d_band_width_Ru_eV": d_width_Ru,
        "d_band_center_Cu_eV": d_center_Cu,
        "d_band_width_Cu_eV": d_width_Cu,
    }

    return record


def find_screening_confs(root_dir):
    """
    Locates all calculation folders matching:
    Intermediate/dock/screening/conf_* or Intermediate/dock/screening_final/conf_*
    Prioritizes 'screening_final' if both exist for an intermediate.
    """
    conf_dirs = []
    abs_root = os.path.abspath(root_dir)

    # 1. Check if user is running directly within a screening folder
    local_confs = sorted(glob.glob(os.path.join(abs_root, "conf_*")))
    valid_local = [c for c in local_confs if os.path.isdir(c)]
    if valid_local:
        return valid_local

    # 2. Walk directory tree looking for dock directories
    for dirpath, dirnames, filenames in os.walk(abs_root):
        if os.path.basename(dirpath) == "dock":
            # Check screening_final first, then screening
            final_dir = os.path.join(dirpath, "screening_final")
            scr_dir = os.path.join(dirpath, "screening")
            target_dir = None
            if os.path.isdir(final_dir):
                target_dir = final_dir
            elif os.path.isdir(scr_dir):
                target_dir = scr_dir

            if target_dir:
                confs = sorted(glob.glob(os.path.join(target_dir, "conf_*")))
                for c in confs:
                    if os.path.isdir(c):
                        conf_dirs.append(c)

    return sorted(list(set(conf_dirs)))


def main():
    parser = argparse.ArgumentParser(description="Extract VASP features for intermediate screening.")
    parser.add_argument("--root", type=str, default=".", help="Root directory to search (default: current directory).")
    parser.add_argument("--output", type=str, default="screening_features.csv", help="Output CSV filename.")
    args = parser.parse_args()

    print(f"Scanning directory: {os.path.abspath(args.root)}")
    conf_dirs = find_screening_confs(args.root)
    print(f"Found {len(conf_dirs)} calculation folders.")

    if not conf_dirs:
        print("No calculation directories found.")
        sys.exit(0)

    records = []
    for i, d in enumerate(conf_dirs, 1):
        print(f"[{i}/{len(conf_dirs)}] Processing: {d}")
        rec = extract_features_from_calc(d)
        if rec is not None:
            records.append(rec)
        else:
            print(f"  -> Skipped (missing or incomplete outputs): {d}")

    if records:
        df = pd.DataFrame(records)
        df.to_csv(args.output, index=False)
        print(f"\nSuccessfully extracted {len(records)} records to: {os.path.abspath(args.output)}")
    else:
        print("No valid calculation records extracted.")

if __name__ == "__main__":
    main()
