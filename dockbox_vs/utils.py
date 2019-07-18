import os
import subprocess

def get_number_of_compounds(file_l):
    suffix, ext = os.path.splitext(file_l)
    if ext == '.mol2':
        # count number of lines with @<TRIPOS>ATOM
        nligs = subprocess.check_output('fgrep -c "@<TRIPOS>ATOM" %s'%file_l, shell=True)
    else:
        raise IOError("Extension not recognized for ligand file!")

    nligs = int(nligs)
    return nligs
