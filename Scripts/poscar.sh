#!/bin/bash

# Create a directory to keep things organized (optional)
mkdir -p collected_poscars

# Loop through all directories matching 'conf_*'
for dir in conf_*; do
    # Check if it's a directory and contains a POSCAR
    if [ -d "$dir" ] && [ -f "$dir/POSCAR" ]; then
        # Extract the number from the folder name (e.g., 0 from conf_0)
        conf_num=${dir#conf_}
        
        # Copy and rename
        cp "$dir/POSCAR" "collected_poscars/POSCAR_$conf_num"
        
        echo "Copied $dir/POSCAR to collected_poscars/POSCAR_$conf_num"
    fi
done

echo "Done! All POSCARs are in the 'collected_poscars' folder."
