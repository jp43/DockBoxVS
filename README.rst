**********
DockBox VS
**********

A plugin to perform VS with DockBoX (DBX)

prepare_vs
**********

prepare_vs creates a directory for vs with all the scripts needed to be run 

Levels of screening with prepare_vs (-level option):

* **Level 0 or 1**: one directory per (ligand, target) pair is created as ligID/targetID. Note that if the number of ligands
is too large (> 1000), subdirectories will be created. For example if this number is 300k, ligand with ligID:lig230134 will
be docked on target target002 in the directory vs/lig230000-lig231000/lig230134/target002 (assuming -w option is not modified).

**N.B.**: The maximum of folders per subdirectory is 1000 while no more than one layer of directories is allowed. Thus the maximum
of ligands allowed for levels 0 or 1 is 1M = 1000*1000. Level 0 and 1 do not remove files subsequent to docking. Hence, all the
binding poses as well as their scores and rescores (in case of rescoring) will be kept, consistently with the options provided in
the DBX configuration file.

Further information on levels 0 and 1:

* **Level 0**: recommended when having a small number (ligand, target) pairs (nligands*ntargets)<1000. One job per (ligand, target) 
pair is created (slurm or SGE script depending on the scheduler). Level 0 is also ideal when having many targets but few compounds
as level 1 systematically docks on all targets in a single job and thus will take more time.

* **Level 1**: recommended when having an intermediate number of (ligand, target) pairs, typically (nligands*ntargets)<1M. Multiple
ligands will be docked (on all targets) in a single job (slurm or SGE script depending on the scheduler). The number of ligands
docked per job can be controlled with the -nligands-per-job option. 

* **Level 2**: recommended for a wide VS screening (nligands*ntargets)>1M. One directory per job will be created inside the VS folder.
Each folder contains a single script (slurm, sge) which will run docking sequentially on multiple ligands (on all targets). All docking
simulations will be run in the same folder. Only the score of the best pose (no rescoring) will be saved


:: 

    usage: prepare_vs [-h] [-l FILE] [-r FILE] [-f FILE] [-level INT]
                  [-nligands-per-job INT] [-s FILE] [-w DIRECTORY NAME]
                  [-slurm OPTIONS | -sge OPTIONS]

    Build directories and config files for Virtual Screening (4th stage)

    optional arguments:
      -h, --help            show this help message and exit
      -l FILE               ligand file: .csv (default: compounds.csv)
      -r FILE               target file(s): .csv (default: targets.csv)
      -f FILE               config file: .ini
      -level INT            Level of screening considered. 0: few compounds; 1:
                            intermediate number of compounds, 2: large number of
                            compounds
      -nligands-per-job INT
                            Number of ligands to be run for every submitted job
                            (when level of VS is 1 or 2))
      -s FILE               csvfile with binding sites: .csv (default: sites.csv)
      -w DIRECTORY NAME     name of directory created for virtual screening
      -slurm OPTIONS        Options for Slurm Workload Manager
      -sge OPTIONS          Options for Sun Grid Engine


* **Examples**

* Excluding specific nodes with sge:

::
    prepare_vs -l compounds.csv -level 1 -sge N,vs,l,'h=!sl390lin22' -f config.ini -nligands-per-job 400


