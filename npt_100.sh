#!/bin/bash
#SBATCH -J production_100_test            # job name
#SBATCH -o prod_100_test.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested 6 288
#SBATCH -n 128       # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 36:00:00         # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

set -x                       # Echo commands, use "set echo" with csh
source ~/.bashrc
mkdir -p test/100ns

cp inputfiles/md_100.mdp test/100ns

cd test/production

cp topol_dih2.top md_0_5.gro md_0_5.cpt index.ndx ../100ns/.

cd ../100ns/

if [ ! -f "md_5_100.gro" ]; then

ibrun -np 1 gmx_mpi grompp -f md_100.mdp -c md_0_5.gro -t md_0_5.cpt -p topol_dih2.top -n index.ndx -o md_5_100.tpr

ibrun gmx_mpi mdrun -deffnm md_5_100

ibrun -np 1 gmx_mpi convert-tpr -f md_5_100.tpr -extend 100000 -o md_5_100.tpr -n index.ndx

ibrun gmx_mpi mdrun -deffnm md_5_100 -append -cpi md_5_100.cpt


else

ibrun -np 1 gmx_mpi convert-tpr -f md_5_100.tpr -extend 100000 -o md_5_100.tpr -n index.ndx

ibrun gmx_mpi mdrun -deffnm md_5_100 -append -cpi md_5_100.cpt


fi

