import sys
import setuptools

def exit_with_error(head, body=''):
    _print_admonition('error', head, body)
    sys.exit(1)

# check Python version
if not (sys.version_info[0] == 2 and sys.version_info[1] >= 6):
    exit_with_error("You need Python 2.6.x or Python 2.7.x to install the DockBox package!")

setuptools.setup(name='dockbox-vs',
    version='1.1',
    scripts=['bin/prepare_compounds', 'bin/prepare_targets', 'bin/prepare_sites', 'bin/prepare_vs'],
    install_requires=['dockbox'],
    description='Plug-in for DockBox package',
    long_description=open('README.rst').read(),
)
