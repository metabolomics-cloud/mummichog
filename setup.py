#!/usr/bin/env python

from setuptools import setup, find_packages

with open("asari/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])
    
with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = f.read()
    
    
setup(
  name='mummichog',
  version=__version__,

  author='Shuzhao Li, Francisco Castellanos, Andrei Todor',
  author_email='shuzhao.li@gmail.com',
  description='Pathway and network analysis for metabolomics data',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/shuzhao-li/mummichog',
  license='BSD 3-Clause',
  keywords='metabolomics analysis bioinformatics mass spectrometry systems biology',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=find_packages(include=['mummichog', 'mummichog.*']),
  data_files=[ ('mummichog/resources', ['mummichog/resources/plotly-graphs.js', 'mummichog/resources/plotly-latest.min.js']),
               ('mummichog/tests', ['mummichog/tests/testdata0710.txt']) ],
  include_package_data=True,
  zip_safe=True,
  entry_points = {
        'console_scripts': ['mummichog=mummichog.command_line:main'],
    },

  python_requires='>=3.4',
  install_requires=requirements.splitlines(),

)
