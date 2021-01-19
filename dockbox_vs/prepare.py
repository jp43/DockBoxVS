import os
import sys
import shutil
import argparse
import subprocess
import ConfigParser
from time import time

from glob import glob
import pandas as pd

import queuing
import utils

def extract_content_from_file(filename):
    content = []
    with open(filename, 'r') as ff:
        for line in ff:
            content.append(line)
    return content

def write_content_to_file(filename, content):
    with open(filename, 'w') as ff:
        for line in content:
            ff.write(line)

def build_config_files_for_vs(model, targetids, csvfile):
    """Build config files from model and csvfile file"""

    config_file_basename = os.path.basename(model)
    config_suff, config_ext = os.path.splitext(config_file_basename)

    config_dir = config_suff
    shutil.rmtree(config_dir, ignore_errors=True)
    os.mkdir(config_dir)

    new_config_files = []
    content = extract_content_from_file(model)
    content_no_site = remove_old_site_information(content)

    # add site information
    df_sites = pd.read_csv(csvfile)
    for idx, targetid in enumerate(targetids):
        new_config_file = config_dir + '/' + config_suff + '_%i'%(idx+1) + config_ext
        new_content = add_new_site_information(content_no_site, targetid, df_sites)
        write_content_to_file(new_config_file, new_content)
        new_config_files.append(new_config_file)

    return new_config_files

def remove_old_site_information(config_file_content):
    """Remove section 'SITE' and option site in DOCKING section of config file if exists"""
    isdock = False
    sitesection = False
    docksection = False
    new_content = []

    for line in config_file_content:
        # check if still in section SITE*
        if line.startswith('[SITE'):
            sitesection = True
        if sitesection and line.startswith('[') and not line.startswith('[SITE'): # new section has been reached
            sitesection = False
        # check if still in section DOCKING
        if line.startswith('[DOCKING]'):
            docksection = True
            isdock = True
        if docksection and line.startswith('[') and not line.startswith('[DOCKING]'): # new section has been reached
            docksection = False
        # check if option line in section DOCKING
        if line.strip().startswith('site') and docksection:
            siteline = True
        else:
            siteline = False
        if not sitesection and not siteline:
            new_content.append(line)
    return new_content

def add_new_site_information(config_file_content, targetid, df_sites):
    """Add section 'SITE' and site in DOCKING section according to targetid and site info"""

    rows = df_sites[df_sites['targetID']==targetid]
    if rows.empty:
        sys.exit("No information regarding the site of target %s in .csv file"%targetid)
    nsites = len(rows)

    new_content = []
    if nsites == 1:
        # add new sections 'SITE' and option site
        for line in config_file_content:
            new_content.append(line)
        for idx, row in rows.iterrows():
            new_content.append("\n[SITE]\n")
            new_content.append("center = %s\n"%(row['center']))
            new_content.append("boxsize = %s\n"%(row['size']))
    elif nsites > 1:
        # add new sections 'SITE' and option site
        for line in config_file_content:
            new_content.append(line)
            if line.startswith('[DOCKING]'):
                new_content.append('site = ' + ', '.join(['site%s'%int(row[1]['site']) for row in rows.iterrows()])+'\n')
        for idx, row in rows.iterrows():
            new_content.append('\n[SITE%s]\n'%str(int(row['site'])))
            new_content.append("center = %s\n"%(row['center']))
            new_content.append("boxsize = %s\n"%(row['size']))

    new_content = add_target_information(new_content, targetid)
    return new_content

def add_target_information(content, targetid):
    new_content = []
    for line in content:
        newline = line.replace("$targetid", targetid)
        newline = newline.replace("${targetid}", targetid)
        newline = newline.replace("${target_id}", targetid)
        newline = newline.replace("$target_id", targetid)
        new_content.append(newline)
    return new_content

def check_config_file(model, level):
    suff, ext = os.path.splitext(model)
    if ext != '.ini':
        sys.exit("Config file must be in .ini format!")

    if level == 2:
        check_config_file_lvl2(model)

def check_config_file_lvl2(model):
    config = ConfigParser.SafeConfigParser()
    config.read(model)

    if config.has_option('DOCKING', 'rescoring'):
        yesno = config.get('DOCKING', 'rescoring').lower()
        if yesno == 'yes':
            raise sys.exit("Rescoring not allowed with level 2 of VS!")

    programs = config.get('DOCKING', 'program').split(',')
    if len(programs) != 1:
        sys.exit("Only one docking program with level 2 of VS!")

def is_rescoring(config_file):
    # check config file
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    is_rescoring = False

    # check if rescoring was performed
    if config.has_option('DOCKING', 'rescoring'):
        is_rescoring_config = config.get('DOCKING', 'rescoring').lower()
        if is_rescoring_config in ['yes', 'true']:
            programs = config.get('RESCORING', 'program').lower()
            is_rescoring = True
    return is_rescoring

def get_programs(config_file, is_rescoring=False, nprograms=None):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    section = is_rescoring*'RESCORING'+(not is_rescoring)*'DOCKING'
    programs = config.get(section, 'program').split(',')
    if nprograms is not None and len(programs) != nprograms:
        raise ValueError("Number of programs should be equal to %i!"%nprograms)
    return programs

def get_sites(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    if config.has_option('DOCKING', 'site'):
        sites = map(str.strip, config.get('DOCKING', 'site').split(','))
    else:
        sites = [""]
    return sites

def get_targets_info(csvfile):
    """get information about targets"""
    df_targets = pd.read_csv(csvfile)

    input_files_r = map(lambda x: os.path.abspath(x), list(df_targets.pdbfile))
    targetids = list(df_targets['targetID'])
    return input_files_r, targetids

def get_ligands_iterator(csvfile, chunksize=10000):
    return pd.read_csv(csvfile, iterator=True, chunksize=chunksize)

def get_number_of_ligands(csvfile):
    nligands = int(subprocess.check_output("wc -l %s"%csvfile, shell=True).split()[0])-1
    print "Number of compounds found: %s"%nligands
    return nligands

def get_last_ligid(csvfile):
    ligid = subprocess.check_output("tail -n 1 %s"%csvfile, shell=True).split(',')[0]
    if not ligid.startswith("lig"):
       sys.exit("First column of csvfile %s does not correspond to ligand ID!")
    return ligid

def create_directory(dirname, overwrite=False):
    if overwrite:
        shutil.rmtree(dirname, ignore_errors=True)
    os.mkdir(dirname)

class PrepareVS(object):

    def __init__(self):
        self.nfolders_per_layer = 1000 # number of folders per layer of subdirectory (level 0 or 1)

    def initalize_lvl1(self, nligands):
        nlayers = 0 # initial number of layers
        nfolders_layer = nligands

        # check how many folder layers are needed
        while nfolders_layer/self.nfolders_per_layer > 1:
            nlayers += 1
            nfolders_layer_new = nfolders_layer/self.nfolders_per_layer

            if nfolders_layer%self.nfolders_per_layer > 0:
                nfolders_layer_new += 1
            nfolders_layer = nfolders_layer_new
        return nlayers, ""

    def initalize_lvl2(self, nligands, nligands_per_job):
        nlayers = 0
        njobs = nligands/nligands_per_job
        if nligands%nligands_per_job != 0:
            njobs += 1
        return njobs, 0, ""

    def write_job_scripts(self, vsdir, iterator_l, nligands, input_files_r, targetids, config_files, level, nligands_per_job=None):
        """Write all the scripts needed"""
        # initialize all variables
        jobidx = 0
        ligand_jdx_abs = 0
        if level in [0, 1]:
            nlayers, workdirs = self.initalize_lvl1(nligands)

        elif level == 2:
            njobs, nlayers, script = self.initalize_lvl2(nligands, nligands_per_job)
            program = get_programs(config_files[0], nprograms=1)[0]
            sites = get_sites(config_files[0])

            basename = program
            if len(sites) > 1: basename += '.'

            extract_lines_1 = ""
            extract_lines_2 = ""
            for idx in range(len(sites)):
                dockdir = basename+sites[idx]
                jdx = idx + 1
                extract_lines_1 += "\nscore_%i="%jdx
                extract_lines_2 += "\n  if [[ -f %(dockdir)s/score.out && -f %(dockdir)s/pose-1.mol2 ]]; then\n    \
score_%(jdx)s=`cat %(dockdir)s/score.out`;\n    cat %(dockdir)s/pose-1.mol2 >> poses.mol2;\n  fi"%locals()

        # create directory which contains scripts to be submitted
        submit_dir = 'to_submit_' + vsdir
        create_directory(submit_dir, overwrite=True)

        # iterate over all compounds
        for chunk in iterator_l:
            input_files_l = map(lambda x: os.path.abspath(x), list(chunk.file_origin))

            input_files_l_indices = list(chunk['index_file'])
            ligand_ids = list(chunk['ligID'])
    
            for jdx, ligid in enumerate(ligand_ids):
                ligand_jdx_abs += 1

                # if VS level is 0 or 1, create a subdirectory per (ligand, target) pair
                if level in [0, 1] or level is None:
                    subdir = vsdir
                    for layer_idx in range(nlayers, 0, -1):
                        subdir += utils.get_subdir_from_ligid(ligid, nligands, self.nfolders_per_layer, layer_idx=layer_idx)
                    for idx, targetid in enumerate(targetids):
                        workdir = subdir + '/' + ligid + '/' + targetid
                        if idx == 0:
                            input_file_l_rel = os.path.abspath(input_files_l[jdx])
    
                        if ligand_jdx_abs == 1 and idx == 0:
                            input_files_r_rel = []
                            for file_r in input_files_r:
                                input_files_r_rel.append(os.path.abspath(file_r))
    
                            config_files_rel = []
                            for file_c in config_files:
                                config_files_rel.append(os.path.abspath(file_c))

                        os.makedirs(workdir) # make directory
                        script = ""
                        if level == 0:
                            filename = workdir + "/run." + self.scheduler
                            script += queuing.make_header(self.scheduler_options, self.scheduler)
                        else:
                            filename = workdir + "/run.sh"
                            script += '#!/bin/bash'
                        script += "\n\nsplit_mol2files -i %s -o ligand.mol2 -ix %s"%(input_file_l_rel, input_files_l_indices[jdx])
                        script += "\ncp %s config.ini"%config_files_rel[idx] # copy to avoid memory problems
                        script += "\ncp %s target.pdb"%input_files_r_rel[idx] # copy to avoid memory problems
                        script += "\nrundbx -f config.ini -l ligand.mol2 -r target.pdb"
                        script += "\nrm -rf config.ini ligand.mol2 target.pdb"
                    
                        with open(filename, 'w') as ff:
                            ff.write(script)

                    if level == 1:
                        # update workdirs for scripts
                        workdirs += ' ' + subdir + '/' + ligid + '/target*'
                        # check if job script should be written
                        if ligand_jdx_abs%nligands_per_job == 0 or ligand_jdx_abs == nligands:
                            jobidx += 1
                            script = queuing.make_header(self.scheduler_options, self.scheduler, jobname="vs%s"%jobidx, output="/dev/null", error="/dev/null")
                            script += """\ndirs=`echo %(workdirs)s`
curdir=`pwd`
for dir in $dirs; do
  cd $dir
  bash run.sh
  cd $curdir
done\n"""%locals()
                            with open(submit_dir+'/run_vs_%i.'%jobidx+self.scheduler, 'w') as ff:
                                ff.write(script)
                            workdirs = ""
 
                elif level == 2:
                    for idx, targetid in enumerate(targetids):
                        if idx == 0:
                            input_file_l_rel = os.path.abspath(input_files_l[jdx])
  
                        if ligand_jdx_abs == 1 and idx == 0:
                            input_files_r_rel = []
                            for file_r in input_files_r:
                                input_files_r_rel.append(os.path.abspath(file_r))
  
                            config_files_rel = []
                            for file_c in config_files:
                                config_files_rel.append(os.path.abspath(file_c))
  
                        mol2cmd = "split_mol2files -i %s -o ligand.mol2 -ix %s"%(input_file_l_rel, input_files_l_indices[jdx])
                        cf = config_files_rel[idx]
                        rf = input_files_r_rel[idx]

                        extract_lines_3 = ""
                        for kdx in range(len(sites)):
                            ldx = kdx + 1
                            extract_lines_3 += "\necho \"%(ligid)s,%(targetid)s,$score_%(ldx)i,%(ldx)i\" >> scores.dat"%locals()
                        script += """\n\n%(mol2cmd)s%(extract_lines_1)s
if [ -f ligand.mol2 ]; then
  cp %(cf)s config.ini
  cp %(rf)s target.pdb
  rundbx -f config.ini -l ligand.mol2 -r target.pdb%(extract_lines_2)s
fi%(extract_lines_3)s
rm -rf config.ini ligand.mol2 target.pdb %(program)s* poses\n"""%locals()
                    # check if job script should be written
                    if ligand_jdx_abs%nligands_per_job == 0 or ligand_jdx_abs == nligands:
                        jobidx += 1
                        jobid = '0'*(len(str(njobs))-len((str(jobidx)))) + str(jobidx)
  
                        jobdir = vsdir + "/job%s"%jobid
                        os.mkdir(jobdir)
  
                        script = '\n\necho "ligID,targetID,%(program)s,site" > scores.dat'%locals() + script
                        script = queuing.make_header(self.scheduler_options, self.scheduler, jobname="vs%s"%jobidx, \
                                output=self.scheduler+"-vs-%s.out"%jobidx, error=self.scheduler+"-vs-%s.err"%jobidx) + script
                        with open(jobdir+'/run.'+self.scheduler, 'w') as ff:
                            ff.write(script)
                        script = ""

        self.write_submit_all_script(submit_dir, vsdir, level, nlayers)

    def write_submit_all_script(self, submit_dir, vsdir, level, nlayers):
        script = submit_dir + '/submit_all.sh'
        scheduler = self.scheduler
        exe = self.exe
        if level == 0:
            dirs = vsdir
            for idx in range(nlayers):
                dirs += '/lig*'
            dirs += '/lig*/target*'
            line_iterate = "for dir in %s; do"%dirs
            # write script to submit all jobs
            script_all = """#!/bin/bash
curdir=`pwd`
%(line_iterate)s
  cd $dir
  %(exe)s run.%(scheduler)s
  cd $curdir
done\n"""%locals()
            with open(script, 'w') as ff:
                ff.write(script_all)
        
        elif level == 1:
            with open(script, 'w') as ff:
                ff.write("""#!/bin/bash
for file in %(submit_dir)s/run_vs_*.%(scheduler)s; do
  %(exe)s $file
done"""%locals())

        elif level == 2:
            with open(script, 'w') as ff:
                ff.write("""#!/bin/bash
curdir=`pwd`
for dir in %(vsdir)s/job*; do
  cd $dir
  %(exe)s run.%(scheduler)s
  cd $curdir
done"""%locals())

    def write_logfile(self, filename, csvfile_l, csvfile_r, csvfile_s, config_file, level):

        with open(filename, 'w') as logf:
            script = """# Virtual screening logfile

# configuration file: %(config_file)s
# ligand file: %(csvfile_l)s
# receptor file: %(csvfile_r)s
# sites file: %(csvfile_s)s
# screening level: %(level)s"""%locals()
            logf.write(script)

    def prepare_scripts(self, workdir, csvfile_l, csvfile_r, csvfile_s, config_file, level, nligands_per_job=None):
        """Prepare all the scripts needed for VS"""
        # extract information about targets (file locations, IDs)
        input_files_r, targetids = get_targets_info(csvfile_r)
    
        # prepare config files
        check_config_file(config_file, level) # check provided model of config file
        config_files = build_config_files_for_vs(config_file, targetids, csvfile_s)
    
        # get csvfile iterator for ligands (needed when dealing with many compounds)
        nligands = int(get_last_ligid(csvfile_l)[3:])

        iterator_l = get_ligands_iterator(csvfile_l, chunksize=10000)
        if level is not None: 
            # write scripts
            create_directory(workdir, overwrite=True)
            self.write_job_scripts(workdir, iterator_l, nligands, input_files_r, targetids, config_files, level, nligands_per_job=nligands_per_job)

    def set_scheduler(self, args, level):
        self.scheduler = None
        for sch in queuing.known_schedulers:
            args_dict = vars(args)
            if args_dict[sch+'_options'] is not None:
                self.scheduler = sch
                self.exe = queuing.exes[sch]
                self.scheduler_options = queuing.check_scheduler_options(args_dict[sch+'_options'], self.scheduler)

        if self.scheduler is None:
            self.exe = None
            self.scheduler_options = None
            if level is not None:
                sys.exit('Info about the scheduler should be provided when level option is specified!')

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="Build directories and config files for Virtual Screening (4th stage)")
        
        parser.add_argument('-l',
            type=str,
            dest='csvfile_l',
            metavar='FILE',
            default='compounds.csv',
            help='ligand file: .csv (default: compounds.csv)')
        
        parser.add_argument('-r',
            type=str,
            dest='csvfile_r',
            metavar='FILE',
            default='targets.csv',
            help = 'target file(s): .csv (default: targets.csv)')
        
        parser.add_argument('-f',
            type=str,
            dest='config_file',
            metavar='FILE',
            default='config.ini',
            help='config file: .ini')
        
        parser.add_argument('-level',
            dest='vs_level',
            type=int,
            metavar='INT',
            choices=range(3),
            help='level of screening considered. 0: few compounds; 1: intermediate number of compounds, 2: large number of compounds')
        
        parser.add_argument('-nligands-per-job',
            dest='nligands_per_job',
            type=int,
            metavar='INT',
            help='number of ligands to be run for every submitted job (should be used with VS level 1 or 2))')

        parser.add_argument('-s',
            type=str,
            dest='csvfile_s',
            metavar='FILE',
            default='sites.csv',
            help = 'sites file(s): .csv (default: sites.csv)')

        parser.add_argument('-w',
            dest='vsdir',
            type=str,
            default='vs',
            metavar='DIRECTORY NAME',
            help='name of directory created for virtual screening')
        
        group = parser.add_mutually_exclusive_group(required=False)
        for sch in queuing.known_schedulers:
            group.add_argument('-%s'%sch,
                dest='%s_options'%sch,
                type=str,
                metavar='OPTIONS',
                help='Options for %s'%queuing.known_schedulers[sch])
        return parser

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()

        if args.vs_level in [1, 2] and not args.nligands_per_job:
            sys.exit('-nligands-per-job option (number of ligands per job) is required when VS level is 1 or 2!')

        self.set_scheduler(args, args.vs_level)
        self.prepare_scripts(args.vsdir, args.csvfile_l, args.csvfile_r, args.csvfile_s, args.config_file, args.vs_level, \
nligands_per_job=args.nligands_per_job)

        # write .log file
        self.write_logfile("config_%s.log"%args.vsdir, args.csvfile_l, args.csvfile_r, args.csvfile_s, args.config_file, args.vs_level)

class CheckVS(PrepareVS):

    def create_arg_parser(self):
        #TODO: add the possibilty to prepare scripts for torque and ll
        parser = argparse.ArgumentParser(description="Check VS and possibly rebuild folders for unfinished ligands!")

        parser.add_argument('--extract',
            dest='extract',
            action='store_true',
            default=False,
            help='Extract results!')
       
        parser.add_argument('-log',
            type=str,
            dest='logfile',
            metavar='FILE',
            default='config_vs.log',
            help='logfile containing info about screening (default: config_vs.log)')

        parser.add_argument('-ntops',
            dest='ntops',
            type=int,
            metavar='INT',
            default=1000,
            help='Number of top hits saved (default: 1000)')

        parser.add_argument('--resume',
            dest='resume',
            action='store_true',
            default=False,
            help='Prepare scripts for unfinished jobs (default: False)')

        parser.add_argument('-vsdir',
            type=str,
            dest='vsdir',
            metavar='FILE',
            default='vs',
            help='VS directory where to extract results from (default: vs)')

        # options for resubmission of unfinished compounds
        parser.add_argument('-nligands-per-job',
            dest='nligands_per_job',
            type=int,
            metavar='INT',
            help='Number of ligands to be run for every submitted job (used for VS levels 1 and 2))')
        
        group = parser.add_mutually_exclusive_group(required=False)
        for sch in queuing.known_schedulers:
            group.add_argument('-%s'%sch,
                dest='%s_options'%sch,
                type=str,
                metavar='OPTIONS',
                default=None,
                help='Options for %s'%queuing.known_schedulers[sch])

        return parser

    def extract_from_logfile(self, filename, ft_name, ftype=str):
        """Extract parameter from log file"""
        with open(filename, 'r') as logf:
            for line in logf:
                line_st = line.strip()
                if line_st.startswith("#"):
                    line_st_sp = line_st[1:].strip().split(":")
                    if len(line_st_sp) == 2:
                        ft = line_st_sp[0].strip()
                        value = line_st_sp[1].strip()
                        if ft == ft_name:
                            return ftype(value)
        sys.exit("Parameter matching description '%s' not found"%ft_name)

    def check_results(self, vsdir, csvfile_l, csvfile_r, config_file, level, extract=False):
        """Check results from previous screening simulations"""

        if extract:
            results_dir = 'results'
            shutil.rmtree(results_dir, ignore_errors=True)
            os.mkdir(results_dir)
        
        # check how many targets per ligand
        ntargets_per_ligands = len(pd.read_csv(csvfile_r))

        if level == 1:
            # extract info about ligands and targets from csvfile
            nligands = len(pd.read_csv(csvfile_l))
            targetID = pd.read_csv(csvfile_r).iloc[-1]['targetID']

            # directory that should be checked
            to_be_checked = 'rescoring/dsx.score'

            nlayers = 0
            # count up how many layers
            subdirs = glob(vsdir+'/lig*')
            while len(subdirs) != 0:
                nlayers += 1
                subdirs = glob(subdirs[0]+'/lig*')
            nlayers = nlayers - 1

            # get subdir names
            subdirs = subprocess.check_output('echo '+vsdir+nlayers*'/lig*', shell=True).split()
            #subdirs = subdirs[:2]

            nligands_done = 0
            dirs_undone = []

            for idx, subdir in enumerate(subdirs):
                dirs_done =  subprocess.check_output('echo %s/lig*/%s/%s '%(subdir, targetID, to_be_checked), shell=True).split()
                dirs = subprocess.check_output('echo %s/lig*'%subdir, shell=True).split()

                dirs_done = ['/'.join(dir.split('/')[:-len(to_be_checked.split("/"))-1]) for dir in dirs_done]
                dirs_undone.extend(list(set(dirs) - set(dirs_done)))
                nligands_done += len(dirs_done)

            # Print info on ligands
            print "Percentage of ligands done: %.2f%% (of %i)"%(nligands_done*100./nligands, nligands)
            print "%i compounds remaining!"%(nligands-nligands_done)
            nligands_undone = len(dirs_undone)
            return sorted(dirs_undone)

        elif level == 2:
            programs = get_programs(config_file, is_rescoring=False)

            results_csv = "%s/scores.csv"%results_dir
            score_files = ['%s/job*/scores.dat'%vsdir for vsdir in vsdirs]
            score_files = ' '.join(score_files)

            # create array with files
            merge_csvfiles_lines = """files=%(score_files)s
first=`echo $files | head -n1 | awk '{print $1;}'`
head -n 1 $first > %(results_csv)s
for f in $files; do
  sed 1d $f >> %(results_csv)s
done
head -1 %(results_csv)s > tmp.csv
sed 1d %(results_csv)s | sort >> tmp.csv
mv tmp.csv %(results_csv)s\n"""%locals()
            subprocess.check_call(merge_csvfiles_lines, shell=True)
            # extract info about ligands from csvfile
            info_l = pd.read_csv(csvfile_l)
            nligands = len(info_l)

            # extract results on ligands
            df = pd.read_csv(results_csv)
            nligands_done = len(df)

            # Print info on ligands
            print "Percentage of ligands done: %.2f%% (of %i)"%(nligands_done*100./nligands/ntargets_per_ligands, nligands)
            print "%i compounds remaining!"%(nligands-nligands_done/ntargets_per_ligands)

            ntops = min(ntops, nligands_done)
            df_groupby = df.groupby('ligID')

            for prgm in programs:
                df_best_score = df_groupby[prgm].idxmin()
                df_best_score_na = df_best_score.isnull()

                df_failed = pd.DataFrame(df_best_score_na[df_best_score_na].index)

                df_no_failed = df.loc[df_best_score[~df_best_score_na]]
                df_no_failed = df_no_failed.reset_index(drop=True)

                # save done ligands
                merged_no_failed = pd.merge(info_l, df_no_failed, on='ligID', how='right')
                merged_no_failed.to_csv("%s/ligands_%s.csv"%(results_dir,prgm), index=False)
                merged_no_failed_smallest = merged_no_failed.nsmallest(ntops, prgm)

                # get all the files with poses 
                filenames = []
                for vsdir in vsdirs:
                    filenames.extend(sorted(list(glob('%s/job*/poses.mol2'%vsdir))))

                poses = {}
                top_ligands = list(merged_no_failed_smallest['name'].values)
                # collect poses if they exist
                for mol2file in filenames:
                    with open(mol2file, 'r') as mol2f:
                        for line in mol2f:
                            if line.startswith("@<TRIPOS>MOLECULE"):
                                name = mol2f.next().strip()
                                if name in top_ligands:
                                    poses[name] = []
                                    poses[name].append(line)
                                    poses[name].append(name+'\n')
                                    keep_molecule = True
                                else:
                                    keep_molecule = False
                            elif keep_molecule:
                                poses[name].append(line)

                merged_no_failed_smallest.to_csv("%s/scores_%s_top_%i.csv"%(results_dir, prgm, ntops), index=False)
                # save failed compounds
                merged_failed = pd.merge(info_l, df_failed, on='ligID', how='right')
                merged_failed.to_csv("%s/ligands_failed_%s.csv"%(results_dir, prgm), index=False)

            merged = pd.merge(info_l, df, on='ligID', how='left', indicator=True)
            undone_ligands = info_l[merged['_merge']=='left_only']

            undone_csv = "%s/ligands_undone.csv"%results_dir
            undone_ligands.to_csv(undone_csv, index=False)
            return undone_csv

    def prepare_scripts_undone_jobs(self, dirs, vsdir, level, nligands_per_job):

        if level == 1:
            jobidx = 0
            nligands = len(dirs)

            # create directory which contains scripts to be submitted
            submit_dir = 'to_submit_' + vsdir
            create_directory(submit_dir, overwrite=True)

            dirs_undone_str = []
            for idx, dir in enumerate(dirs):
                dirs_undone_str.append(dir+'/target*')
                jdx = idx + 1
                if jdx%nligands_per_job == 0 or jdx == nligands:
                    jobidx += 1
                    workdirs = ' '.join(dirs_undone_str)
                    script = queuing.make_header(self.scheduler_options, self.scheduler, jobname="vs%s"%jobidx, output="/dev/null", error="/dev/null")
                    script += """\ndirs=`echo %(workdirs)s`
curdir=`pwd`
for dir in $dirs; do
  cd $dir
  bash run.sh
  cd $curdir
done\n"""%locals()
                    with open(submit_dir+'/run_vs_%i.'%jobidx+self.scheduler, 'w') as ff:
                        ff.write(script)
                    dirs_undone_str = []
        elif level == 2:
            raise NotImplementedError

    def run(self):
        parser = self.create_arg_parser()
        args = parser.parse_args()

        logfile = args.logfile
        # load info about receptor and ligand files
        csvfile_l = self.extract_from_logfile(logfile, 'ligand file')
        csvfile_r = self.extract_from_logfile(logfile, 'receptor file')

        # load info about config file
        config_file = self.extract_from_logfile(logfile, 'configuration file')
        csvfile_s = self.extract_from_logfile(logfile, 'sites file')

        # VS level
        level = self.extract_from_logfile(logfile, 'screening level', ftype=int) 
        dirs_undone = self.check_results(args.vsdir, csvfile_l, csvfile_r, config_file, level, extract=args.extract)

        if args.resume:
            self.set_scheduler(args, level)
            self.prepare_scripts_undone_jobs(dirs_undone, args.vsdir, level, args.nligands_per_job)
