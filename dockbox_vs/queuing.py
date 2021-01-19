import sys
import shlex

known_schedulers = {'sge': 'Sun Grid Engine', 'slurm': 'Slurm Workload Manager'}
exes = {'sge': 'qsub', 'slurm': 'sbatch'}

# default and mandatory options for schedulers
default_options = {'sge': {'S': '/bin/bash',
'o': 'sge-$JOB_ID.out',
'e': 'sge-$JOB_ID.err',
'N': 'vs',
'cwd': '',
'V': ''},
'slurm': {'job-name': 'vs',
'ntasks': '1',
'cpus-per-task': '1',
'nodes': '1'}}

#mandatory_options = {'sge': ['q'], 'slurm': ['partition', 'time']}
mandatory_options = {'sge': [], 'slurm': ['partition', 'time']}
equivalent_options = {'slurm': {'p': 'partition', 't': 'time', 'j': 'job-name', 'x': 'exclude'}}

def slurm_to_seconds(string):

    # check if days are provided
    string_s = string.split("-")
    if len(string_s) == 2:
        days = int(string_s[0])
        string_s_s = string_s[1].split(":")
        if len(string_s_s) == 2:
            hours = int(string_s_s[0])
            minutes = int(string_s_s[1])
        elif len(string_s_s) == 1:
            hours = int(string_s_s[0])
            minutes = 0
        else:
            raise Exception("SLURM time format %s not recognized"%string)
        seconds = 0
    elif len(string_s) == 1:
        string_s_s = string_s[0].split(":")
        days = 0
        if len(string_s_s) == 3:
            hours = int(string_s_s[0])
            minutes = int(string_s_s[1])
            seconds = int(string_s_s[2])
        elif len(string_s_s) == 2:
            hours = int(string_s_s[0])
            minutes = int(string_s_s[1])
            seconds = 0
        else:
            raise Exception("SLURM time format %s not recognized"%string)

    time = days*24*3600 + hours*3600 + minutes*60 + seconds
    return time

def check_scheduler_options(string, scheduler):

    options = shlex.shlex(string, posix=True)
    options.whitespace += ','
    options.whitespace_split = True

    options = list(options)
    noptions = len(options)/2
    if len(options) != 2*noptions:
        raise IOError('Number of options provided for th scheduler should be even!')

    options_dict = {}
    for key, value in default_options[scheduler].iteritems():
        options_dict[key] = value

    for idx in range(noptions):
        key = str(options[2*idx])
        value = str(options[2*idx+1])
        if scheduler in equivalent_options and key in equivalent_options[scheduler]:
            key = equivalent_options[scheduler][key]
        options_dict[key] = value

    # check if mandatory options have been set
    for mo in mandatory_options[scheduler]:
        if mo not in options_dict.keys():
            raise ValueError('Option %s mandatory for %s scheduler is not provided'%(mo,scheduler))
    return options_dict

def make_header(options, scheduler, jobname=None, output=None, error=None):

    header = "#!/bin/bash"
    for key, value in options.iteritems():
        if scheduler.lower() == 'sge':
            if jobname is not None and key == 'N':
                value = jobname
            elif output is not None and key == 'o':
                value = output
            elif error is not None and key == 'e':
                value = error
            header += '\n#$ -%s %s'%(key, value)
        elif scheduler.lower() == 'slurm':
            if len(key) == 1:
                if jobname is not None and key == 'j':
                    value = jobname
                header += '\n#SBATCH -%s=%s'%(key, value)
            else:
                if jobname is not None and key == 'job-name':
                    value = jobname
                header += '\n#SBATCH --%s=%s'%(key, value)
    return header
