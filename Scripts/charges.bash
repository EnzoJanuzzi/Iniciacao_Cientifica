#!/bin/bash

# Define the parent directory name
PARENT="Charges"

# Define the 3 sub-directories
SUBDIRS=("Sistema" "NP" "Molec")

# Loop through all site folders (e.g., atop1/, bridge/, etc.)
for site in */; do
    # Remove trailing slash for cleaner echo outputs
    site_clean="${site%/}"

    echo "---------------------------------------------------------"
    echo " Processing: $site_clean"

    # Enter the site directory
    cd "$site" || continue

    # Check if CONTCAR exists inside this specific site folder
    if [ ! -f "CONTCAR" ]; then
        echo " [ERROR] CONTCAR not found in $site_clean! Skipping..."
        cd ..
        continue
    fi

    # Check if Charges folder doesn't exist
    if [ ! -d "$PARENT" ]; then
        mkdir -p "$PARENT"

        # Loop through the array to create subfolders and copy files
        for DIR in "${SUBDIRS[@]}"; do
            TARGET="$PARENT/$DIR"

            # Create the subfolder
            mkdir -p "$TARGET"

            # Copy CONTCAR and rename it to POSCAR
            cp CONTCAR "$TARGET/POSCAR"

            # Copy standard VASP files from the site root if they exist
            for FILE in KPOINTS job-vasp.sh job-vasp_final.sh; do
                if [ -f "$FILE" ]; then
                    cp "$FILE" "$TARGET/"
                fi
            done
	    
	    # Copy INCAR charge file
	    cp ~/IC/Intermediarios/Geral/INCAR/INCAR_charges "$TARGET/INCAR"
        done
        echo " [+] Folder '$PARENT' created with subdirectories in $site_clean."
    else
        echo " [!] Folder '$PARENT' already exists in $site_clean (skipping)."
    fi

    # Return to the parent directory at the end of the iteration
    cd ..
done

echo "========================================================="
echo "Setup complete! Folders are ready in each site."
echo "========================================================="
echo "??  Next Manual Steps:"
echo " 1. Edit POSCAR in NP (delete adsorbate atoms & velocities)."
echo " 2. Edit POSCAR in Molec (delete NP atoms & velocities)."
echo " 3. Update the POTCAR files in NP and Molec to match the new POSCARs."
echo " 4. Ensure all INCARs have IBRION = -1, NSW = 0, LAECHG = .TRUE., and LCHARG = .TRUE."
