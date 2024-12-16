#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=scavenge
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=00:04:00
#SBATCH --mem=256G
#SBATCH --nodes=1 
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --mail-user=ht33@rice.edu
#SBATCH --mail-type=ALL   

echo "Loading modules..."
module load GCCcore/13.2.0
module load Python/3.11.5

echo "Copying data to scratch..."
cp -r /projects/comp416/GeneGorman/hard_data /scratch/ht33

echo "Changing directory to code location..."
cd /projects/comp416/GeneGorman

echo "Running main.py..."
./main.py /scratch/ht33/hard_data

echo "Cleaning up..."
rm -r /scratch/ht33/hard_data