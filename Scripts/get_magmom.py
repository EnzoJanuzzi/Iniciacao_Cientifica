#!/usr/bin/env python3
import sys

def extract_magmom(file_path):
    """
    Parses a VASP OUTCAR or text file containing the magnetization (x) table
    and returns a formatted MAGMOM string using the 'tot' column.
    """
    tot_values = []
    inside_table = False

    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Look for the start of the magnetization table
                if "magnetization (x)" in line or "magnetization (x)" in line.lower():
                    inside_table = True
                    tot_values = []  # Reset to get the LAST occurrence in OUTCAR
                    continue

                if inside_table:
                    # Look for the header line separation
                    if "---" in line:
                        continue
                    
                    # Stop reading when reaching the end of the table
                    if "tot" in line or "NOUT" in line or line.strip() == "":
                        if len(tot_values) > 0:
                            inside_table = False
                            continue

                    parts = line.split()
                    # A valid atom line has at least 5 elements: ion_# s p d tot
                    if len(parts) >= 5 and parts[0].isdigit():
                        tot_values.append(parts[-1])

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

    if not tot_values:
        print("Error: Could not find any magnetization table or 'tot' values in the file.")
        sys.exit(1)

    magmom_str = f"MAGMOM = {' '.join(tot_values)}"
    return magmom_str

if __name__ == "__main__":
    # If a file path is provided as argument, use it; otherwise default to 'OUTCAR'
    target_file = sys.argv[1] if len(sys.argv) > 1 else "OUTCAR"
    
    print("\n--- Extracted MAGMOM ---")
    print(extract_magmom(target_file))
    print("------------------------\n")
