"""
Build mirge.
"""
import sys
import os

from setuptools import setup, Extension, find_packages
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext

if sys.version_info < (2, 7):
	sys.stdout.write("At least Python 2.7 is required.\n")
	sys.exit(1)

# Get the long description from the README file
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'README.md')) as f:
	long_description = f.read()

setup(
	name='mirge',
	version='2.0',
	author='Yin Lu',
	author_email='ylu61@jhmi.edu',
	#url='https://github...',
	description='comprehensive analysis of miRNA sequencing data',
	long_description=long_description,
	license='MIT',
	package_dir={'': 'src'},
	packages=find_packages('src', exclude=['.txt']),
	package_data = {'mirge':['models/*.pkl', 'models/*.txt', 'rScripts/*.R']},
	install_requires=['cutadapt==1.11', 'biopython==1.68', 'numpy==1.11.3',
	'scipy==0.17.0', 'matplotlib==2.1.1', 'pandas==0.21.1','scikit-learn==0.18.1',
	'reportlab==3.3.0', 'forgi==0.20',
	],
	entry_points={'console_scripts': ['miRge2.0 = mirge.__main__:main']},
	classifiers=[
		"Development Status :: 1 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
