import os
import subprocess

known_formats = ['.sdf', '.smi', '.mol2']

def get_number_of_compounds(file_l):
    suffix, ext = os.path.splitext(file_l)
    if ext == '.mol2':
        # count number of lines with @<TRIPOS>ATOM
        nligs = subprocess.check_output('fgrep -c "@<TRIPOS>ATOM" %s'%file_l, shell=True)
    else:
        raise IOError("Extension not recognized for ligand file!")

    nligs = int(nligs)
    return nligs

def get_babel_command(file_l, index=None):

    suffix, ext = os.path.splitext(file_l)
    extout = '.mol2'
    for fmt in known_formats:
        if ext == fmt:
            input_format_flag = '-i' + ext[1:]
            output_format_flag = '-o' + extout[1:]
            if fmt in ['.smi', '.sdf']:
                other_flags = '--gen3D'
            else:
                other_flags = ''
    return 'babel %s %s %s ligand.mol2 -f %s -l %s %s 2>/dev/null'%(input_format_flag, file_l, output_format_flag, index, index, other_flags)

