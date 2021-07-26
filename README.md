# MASSA Algorithm
MASSA Algorithm: A tool for separating data sets of molecules into training and test sets. Developed with the objective of preparing data sets for the generation of prediction models in chemoinformatics.

## Instalation
MASSA Algorithm can be installed using pip:
```
pip install MASSA_Algorithm-0.1-py3-none-any.whl
```
To upgrade to the latest version (recommended), also use pip:
```
pip install --upgrade MASSA_Algorithm
```
Alternatively, you can build the latest development version from source:
```
git clone https://github.com/gcverissimo/MASSA_Algorithm.git
cd MASSA_Algorithm
python setup.py install
```
### Requirements
* python: >= 3.7;
* rdkit (RDKit will not be installed automatically with the package. Recommended cross-platform installation via conda. For more information on how to install RDKit, see: https://www.rdkit.org/docs/Install.html);
* numpy;
* pandas;
* matplotlib: >= 3.2;
* scipy: >= 1.6;
* scikit-learn: >= 0.24;
* kmodes:¹ >= 0.10.

## Usage
Once installed, the program can be run directly from the command line:
```
MASSA_Algorithm -i <input_file>.sdf -o <output_file>.sdf
```

A list of optional arguments include:
* Percentage of molecules in training set: ```-p``` or ```--percentage_of_training```
* Number of biological activities for separation: ```-b``` or ```--number_of_biological```
* Name of biological activities for separation: ```-s``` or ```--the_biological_activities```
* Number of principal components in PCA: ```-n``` or ```--number_of_PCs```
* SVD solver parameter for PCA: ```-v``` or ```--svd_solver_for_PCA```
* Extension of image files: ```-t``` or ```--image_type```
* Font size for X-axis of dendrograms: ```-d``` or ```--dendrogram_Xfont_size```
* Font size for X-axis of bar plots: ```-x``` or ```--barplot_Xfont_size```

### Command line help
A full description of the arguments can also be viewed directly from the command line using the command:
```
MASSA_Algorithm -h
```
or
```
MASSA_Algorithm --help
```

## Cite
```
@Misc{veríssimo2021,
    author = {Gabriel Corrêa Veríssimo},
    title = {MASSA Algorithm: Molecular data set sampling for training-test separation},
    howpublished = {\url{https://github.com/gcverissimo/MASSA_Algorithm}},
    year = {2021}
  }
```

## References
[1]: DE VOS, N. J. kmodes categorical clustering library. https://github.com/nicodv/kmodes. 2015-2021.
