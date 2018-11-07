from setuptools import setup, find_packages
import sys
sys.path.append('../atomutil')
sys.path.append('./test')

setup(
    name='atomutil',
    version='0.0.1',
    description='Handle file for POSCAR and xyz',
    # install_requires=['numpy', 'os', 'copy', 'subprocess', 'pkgutil'],
    packages=find_packages(),
    package_data={'atomtools':['*.txt']},
    test_suite='sample_test.suite'
)
