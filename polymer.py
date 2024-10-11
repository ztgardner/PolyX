import os
import subprocess

BOSSDIR= "/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/boss"
CONDAPATH= "/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/miniconda3/bin/activate"
SINGPATH= "/project/cmri235_uksr/shasanka_conda_boss/sla296/singularity/Fast/f.sif"
def make_itp(path_toPDB, name, charge=0):
    name='Base'
    charge = int(charge) or 0
    path_to_pdb = os.path.join(path_toPDB, f"Base.pdb")
    conda_activate = f"source {CONDAPATH} && conda activate ligpg"
    export_bossdir = f" export BOSSdir={BOSSDIR}"
    ligpargen_cmd = f"ligpargen -i {path_to_pdb} -n Base -p {path_toPDB}/Delete -r {name[:3]} -c {charge} -o 0 -cgen CM1A"
    singularity_container = SINGPATH
    obabel_cmd = f"obabel -ipdb {path_to_pdb} -oxyz -O {os.path.join(path_toPDB, f'{name[:3]}.xyz')} "
    cmd = ["singularity", "exec", singularity_container, "bash", "-c",
           f'{conda_activate} && {export_bossdir} && {ligpargen_cmd} && {obabel_cmd}']
    try:
        print("in try")
        print(cmd)
        subprocess.Popen(cmd).wait()

    except subprocess.CalledProcessError:
        print(
            f'Ligpargen was not able to find a parameter, user input files is being used. Please rerun with your own itp and pdb files.')

    move_gro_PDB_and_itp_command= f"mv {os.path.join(path_toPDB, 'Delete',f'Base.gmx.itp {path_toPDB}/{name}.itp')} && rm -rf {path_toPDB}/Delete"
    subprocess.run(move_gro_PDB_and_itp_command,shell=True,check=True)
    return os.path.join(path_toPDB, f'{name}.itp')


