#!/bin/bash

#SBATCH -t 1:00:00                   #Time for the job to run
#SBATCH -J Dihedral_Scan                  #Name of the job
#SBATCH -o Dihedral_Scan_job.o




#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition=CAL48M192_D  #SKY32M192_D
#SBATCH --mail-type ALL
#SBATCH --mail-user nclo224@uky.edu

#SBATCH --account=col_cmri235_uksr

# Source Load Modules
source /home/nclo224/Load_groamcs.sh
module load ccs/conda/intelpython3_full-2020.0
module load ccs/conda/scipy-1.3.0
module load ccs/singularity


#Set up run directory:
mkdir Dihedral_Scan
cd Dihedral_Scan
cp -r ~/Dihedral_Scan_Code/* . #Needs path to DH_scan Code (could also be a git clone?)
cp ../Opt.log . #For Finding the Charges
cp ../Scan.log . #For building everything else


#Prepare base.itp
python Codes/make_base_pdb.py Opt.log Base #Builds the .pdb to feed ligpargen

source /project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/miniconda3/bin/activate
python -c 'import t; t.make_itp("wat_Solvent")'

rm -r __pycache__ bas.xyz Base.pdb t.py


 

    




#Build Needed Dir
mkdir GRO
mkdir GJF
mkdir ITP





for ((i=0; i<3;i++)); do #Running 3 Times
    python dihedral_scan.py 
    
    largest_run_dir=$(ls -d run_* | sort -r -k 6,7 | head -n 1) #Finds the largest run_# dir
    cd "$largest_run_dir"

    for filename in *.gro; do  
        file="${filename%.*}"  
        # Perform editconf
        gmx_mpi editconf -f ${file}.gro -o ${file}_box.gro -c -d 2.0 -bt cubic

        # Perform grompp
        gmx_mpi grompp -f ../em.mdp -c ${file}_box.gro -p ${file}.top -o em_${file}.tpr -maxwarn 5

        # Perform mdrun
        gmx_mpi mdrun -deffnm em_${file} -v &> em_${file}.mdrun -nb cpu
    done
    
    python paramaterization.py 
    cd ../
done


