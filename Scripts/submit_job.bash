dirs=$(ls -d */)
for dir in $dirs; do
    cd $dir
    sbatch job-vasp.sh 
    cd ..
done
