export CONDA_PREFIX=$(conda info --base)  # This will get the path to the base conda environment

for run in {1..1000}; do
    job_name="opt_randtobest1bin_sobol_run-${run}"
    sbatch_file="${job_name}.slurm"

    # Create sbatch file
    echo "#!/bin/bash" > $sbatch_file
    echo "#SBATCH --job-name=f_${run}" >> $sbatch_file
    echo "#SBATCH --output=${job_name}.out" >> $sbatch_file  # Changed to .out to avoid overwriting the sbatch file
    echo "#SBATCH --ntasks=40" >> $sbatch_file
    echo "#SBATCH --time=10:00:00" >> $sbatch_file
    echo "source $CONDA_PREFIX/bin/activate base" >> $sbatch_file  # Activate the base conda environment
    echo "$CONDA_PREFIX/bin/python diff_ev_trial_gen.py randtobest1bin sobol $job_name" >> $sbatch_file

    # Submit job
    sbatch $sbatch_file

    # Optionally, remove the sbatch file after submission
    rm $sbatch_file -f
    
done


