import setuptools

with open('README.md', 'r') as readme:
	long_description = readme.read()

setuptools.setup(
	name='MASSA Algorithm',
	version='0.9.1',
	description='MASSA Algorithm is a Python package to separate data sets of molecules into training and test sets, considering the diversity of structural, physicochemical and biological characteristics of these molecules.',
	long_description=long_description,
	long_description_content_type='text/markdown',
	author='Gabriel Corrêa Veríssimo',
	author_email='gcverissimo@outlook.com',
	url='https://github.com/gcverissimo/MASSA_Algorithm',
	license='AGPLv3',
	packages= setuptools.find_packages(),
	entry_points={'console_scripts':['MASSA_Algorithm = MASSA_Algorithm.MASSA:main']},
	install_requires=['numpy', 'rdkit-pypi', 'pandas', 'matplotlib >= 3.2', 'scipy >= 1.6', 'scikit-learn >= 0.24', 'kmodes >= 0.10'],
	python_requires='>=3.8, <4',
	keywords=['chemoinformatics', 'training', 'test', 'training-test', 'dataset preparation', 'data set preparation', 'rdkit'],
	classifiers=['Programming Language :: Python :: 3',
	'Programming Language :: Python :: 3.8',
	'Programming Language :: Python :: 3.9',
	'Programming Language :: Python :: 3.10',
	'License :: OSI Approved :: GNU Affero General Public License v3',
	'Operating System :: OS Independent',
	'Topic :: Scientific/Engineering :: Chemistry'])