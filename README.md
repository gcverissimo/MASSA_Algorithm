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
* python: >= 3.7
* rdkit (Recommended cross-platform installation via conda. See: https://www.rdkit.org/docs/Install.html)
* numpy
* pandas
* matplotlib: >= 3.2
* scipy: >= 1.6
* scikit-learn: >= 0.24
* kmodes:¹ >= 0.10 

## Usage

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
1: DE VOS, N. J. kmodes categorical clustering library. https://github.com/nicodv/kmodes. 2015-2021.
