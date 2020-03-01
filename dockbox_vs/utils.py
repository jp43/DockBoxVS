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

def prepare_mol2_command(file_l, index=None):

    suffix, ext = os.path.splitext(file_l)
    for fmt in known_formats:
        if ext == fmt:
            input_format_flag = '-i' + ext[1:]
            if fmt in ['.smi', '.sdf']:
                #TODO: check if sdf is 3D in which case, do not generate the 2D structure
                other_flags = '--gen3D' 
            else:
                other_flags = ''

    if index is None:
        cmd = 'babel %s %s -omol2 ligand.mol2 -f %s 2>/dev/null'%(input_format_flag, file_l, other_flags)
    else:
        cmd = 'babel %s %s -omol2 ligand.mol2 -f %s -l %s %s 2>/dev/null'%(input_format_flag, file_l, index, index, other_flags)
    return cmd

def get_subdir_from_ligid(ligid, nligands, nfolders_per_layer, layer_idx=1):
    """get subdirectory name"""

    minbin = (int(ligid[3:])-1)/nfolders_per_layer**layer_idx
    minbin = minbin*nfolders_per_layer**layer_idx + 1

    maxbin = min(minbin+nfolders_per_layer**layer_idx - 1, nligands)

    minbin_str = str(minbin)
    maxbin_str = str(maxbin)

    minbin_str = (len(ligid[3:])-len(minbin_str))*'0'+minbin_str
    maxbin_str = (len(ligid[3:])-len(maxbin_str))*'0'+maxbin_str

    return '/lig' + minbin_str + '-' + maxbin_str

