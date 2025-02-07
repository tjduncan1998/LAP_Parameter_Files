#!/bin/bash
#SBATCH -J test_prep            # job name
#SBATCH -o prep_test.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested 6 288
#SBATCH -n 128        # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 12:00:00         # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

set -x                       # Echo commands, use "set echo" with csh
source ~/.bashrc

# Generate topology for LAP5

wd='test'
mkdir $wd
n_water=30000
n_lipid=389
n_lipid_tot=$(( $n_lipid+1 ))
n_protein=2
x_box=15
y_box=15
z_box=15
pdb_protein='TAZTRPDPP2_channel.pdb'
pdb_lipid='popc_box.gro'
protein_itp='TAZTRPDPP2_channel.itp'


cp inputfiles/topol_dih2.top inputfiles/prot.ndx inputfiles/strong_posre.itp inputfiles/minim* inputfiles/topol_popc.top inputfiles/$pdb_protein inputfiles/$pdb_lipid inputfiles/inflategro.pl $wd
cd $wd

gmx_mpi --version

echo "GENERATING TOPOLOGY FOR LAP5"

#mktop.pl -i $pdb_protein -c charges.txt -o topol_lap5.top -ff opls -conect yes 

#cp topol_lap5.top oldtopol.top

# Replace atom names and add dihedrals for phosphorous

echo "Adding missing dihedrals for LAP5"

#python3 dihed_corr.py

#cp temptopol.top topol_dih2.top

# Include proper forcefields

echo "Adding forcefield paramters"

sed -i 's/oplsaa.ff/oplsaa_mod2_lipid.ff/g' topol_dih2.top

# CHANGE WHICH STRUCTURE YOU ARE REFERENCING

sed "4 a #include \"/work2/07817/tdunca/frontera/gromacs-2020.2/share/top/oplsaa_mod2_lipid.ff/$protein_itp\""  topol_popc.top > topol_dih2.top

sed -i '$d' topol_dih2.top
sed -i '$d' topol_dih2.top
echo "MOL      1" >> topol_dih2.top

#sed '13 a #include "/work2/07817/tdunca/stampede2/gromacs-2020.2/share/top/oplsaa_mod_lipid.ff/tip4p.itp" ' topol_dih2.top > temptopol.top

#sed '13 a #include "/work2/07817/tdunca/stampede2/gromacs-2020.2/share/top/oplsaa_mod_lipid.ff/popc_mod.itp" ' temptopol.top > topol_dih2.top

# Initialize Lipid Bilayer

echo "Adding lipids and water for 5ns equilibration"

#seed=$(( RANDOM % 10000 ))
#
#ibrun -np 1 gmx_mpi insert-molecules -f $pdb_lipid -ci $pdb_lipid -nmol $n_lipid -try 1000 -box $x_box $y_box $z_box -seed $seed -o popc390.gro
#
#ibrun -np 1 gmx_mpi solvate -cp popc390.gro -cs tip4p.gro -o solv.gro -maxsol 30000
#
#echo 'POPC      '$n_lipid_tot >> topol_dopc.top
#echo 'SOL       '$n_water >> topol_dopc.top

# Energy Minimization

ibrun -np 1 gmx_mpi grompp -f minim_lipid.mdp -c $pdb_lipid -p topol_popc.top -o em.tpr

#ibrun -np 1 gmx_mpi mdrun -deffnm em 

echo 0 | ibrun -np 1 gmx_mpi trjconv -s em.tpr -f $pdb_lipid -o popc390_whole.gro -pbc mol -ur compact

REPLY=$(tail -n 1 popc390_whole.gro) 

ibrun -np 1 gmx_mpi editconf -f $pdb_protein -o LAP5_process.gro -c yes

ibrun -np 1 gmx_mpi editconf -f LAP5_process.gro -o LAP5_newbox.gro -c -box $REPLY 

# Minimize LAP5

#ibrun -np 1 gmx_mpi grompp -f minim.mdp -c LAP5_newbox.gro -p topol_dih2.top -o lap5_min.tpr
#
#ibrun gmx_mpi mdrun -v -deffnm lap5_min

#ibrun -np 1 gmx_mpi solvate -cp lap5_min.gro -cs tip4p.gro -o lap5_solv.gro -p topol.top -maxsol 1000

#ibrun -np 1 gmx_mpi grompp -f nvt.mdp -c lap5_solv.gro -p topol_dih2.top -o lap5_nvt.tpr

#ibrun -np 1 gmx_mpi mdrun -v -deffnm lap5_nvt

(echo $REPLY | read $x $y $z)

n_poly=$(head -n 2 LAP5_newbox.gro | tail -n 1)

n_poly=$(( $n_poly+1 ))
echo $n_poly

n_solv=$(head -n 2 popc390_whole.gro | tail -n 1)

echo $n_solv
n_solv=$(( $n_solv+1 ))

n_tot=$(( $n_solv+$n_poly-2 ))

rm system.gro

echo "Lipid + LAP5" >> system.gro
echo '  '$n_tot >> system.gro
tail -n $n_poly LAP5_newbox.gro | head -n $(( $n_poly-1 )) >> system.gro
tail -n $n_solv popc390_whole.gro >> system.gro

sed '5 a #ifdef STRONG_POSRE' topol_dih2.top > temptopol2.top
sed '6 a #include "strong_posre.itp"' temptopol2.top > temptopol3.top
sed '7 a #endif' temptopol3.top > topol_dih2.top

echo ' ' >> topol_dih2.top

tail -n 2 topol_popc.top | head -n 1  >> topol_dih2.top

#echo 0 | ibrun -np 1 gmx_mpi genrestr -f LAP5_newbox.gro -o strong_posre.itp -fc 100000 100000 100000

perl inflategro.pl system.gro 4 POPC 60 system_inflated.gro 5 area.dat

echo Update Topology: `grep ' P1' system_inflated.gro | grep -v UNL | nl | tail -n 1 `

output=$(grep ' P1' system_inflated.gro | grep -v UNL | nl | tail -n 1 | cut -c 4-6)
sed -i "s/200/$output/g" topol_dih2.top

ibrun -np 1 gmx_mpi grompp -f minim_inflategro.mdp -c system_inflated.gro -p topol_dih2.top -r system_inflated.gro -o system_inflated_em.tpr -maxwarn 2

ibrun gmx_mpi mdrun -deffnm system_inflated_em

echo 0 | ibrun -np 1 gmx_mpi trjconv -s system_inflated_em.tpr -f system_inflated_em.gro -o tmp.gro -pbc mol

mv tmp.gro system_shrink1.gro

mkdir system_shrink

cp inflategro.pl strong_posre.itp system_shrink1.gro topol_dih2.top ../inputfiles/water_deletor.pl system_shrink/.

cd system_shrink/


for i in {1..28}

do

perl inflategro.pl system_shrink${i}.gro 0.95 POPC 0 system_shrink$(($i+1)).gro 5 area_shrink1.dat 

done

mv system_shrink$(($i+1)).gro system_shrunk.gro
rm system_shrink*

gmx_mpi editconf -f system_shrunk.gro -o system_shrunk.gro -box `tail -n 1 system_shrunk.gro | awk '{print $1,$2}'` 12 -c yes

ibrun -np 1 gmx_mpi solvate -cp system_shrunk.gro -cs spc216.gro -o system_solv.gro -p topol_dih2.top  

perl water_deletor.pl -in system_solv.gro -out system_solv_fix.gro -ref O3 -middle C27 -nwater 3

new_water=$((`grep SOL system_solv_fix.gro | nl | tail -n 1 | cut -c 2-7` / 3))

sed -i '$ d' topol_dih2.top
echo "SOL $new_water" >> topol_dih2.top

echo 'Update Topology'

cd ../../

mkdir test/equil
cp inputfiles/*.mdp inputfiles/prot.ndx test/equil/ 
cd test/system_shrink/
cp topol_dih2.top system_solv_fix.gro ../strong_posre.itp ../equil/
cd ../equil

ibrun -np 1 gmx_mpi grompp -f minim.mdp -c system_solv_fix.gro -r system_solv_fix.gro -p topol_dih2.top -o em.tpr -maxwarn 2
ibrun gmx_mpi mdrun -deffnm em

#for i in {1..5}; do
#ibrun -np 1 gmx_mpi grompp -f minim_cg.mdp -c em.gro -r em.gro -p topol_dih2.top -o em.tpr -maxwarn 2
#ibrun gmx_mpi mdrun -deffnm em
#done

#echo q | ibrun -np 1 gmx_mpi make_ndx -f em.gro -o index.ndx 
#cat prot.ndx >> index.ndx

ibrun -np 1 gmx_mpi make_ndx -f em.gro -o index.ndx << EOF
2 | 3
q
EOF

ibrun -np 1 gmx_mpi grompp -f nvt_res.mdp -c em.gro -r em.gro -p topol_dih2.top -n index.ndx -o nvtres.tpr -maxwarn 1 
ibrun gmx_mpi mdrun -deffnm nvtres

ibrun -np 1 gmx_mpi grompp -f nvt_eq.mdp -c nvtres.gro -r nvtres.gro -p topol_dih2.top -n index.ndx -o nvt.tpr -maxwarn 1 
ibrun gmx_mpi mdrun -deffnm nvt

ibrun -np 1 gmx_mpi grompp -f npt.mdp -c nvt.gro  -r nvt.gro -t nvt.cpt -p topol_dih2.top -n index.ndx -o npt.tpr -maxwarn 1
ibrun gmx_mpi mdrun -deffnm npt

cd ../../

mkdir -p test/production
cp inputfiles/md.mdp test/production/.
cd test/equil/
cp topol_dih2.top npt.* index.ndx ../production/.
cd ../production

ibrun -np 1 gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol_dih2.top -n index.ndx -o md_0_5.tpr

ibrun gmx_mpi mdrun -deffnm md_0_5
`
