# MASSA Algorithm
MASSA Algorithm: A tool for separating data sets of molecules into training and test sets. Developed with the objective of preparing data sets for the generation of prediction models in cheminformatics.

## Instalation
MASSA Algorithm can be installed using pip:
```
pip install MASSA_Algorithm
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
* python: >= 3.8;
* rdkit;
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
* **Percentage of molecules in training set**: ```-p``` or ```--percentage_of_training```.
    * Percentage of molecules in training set. Must be a number from 0 to 1.
    * Default = 0.8.
* **Number of biological activities for separation**: ```-b``` or ```--number_of_biological```.
    * Number of biological activities that will be used to separate the set into training and test.
    * Default = 1.
* **Name of biological activities for separation**: ```-s``` or ```--the_biological_activities```.
    * Enter a list with the names of biological activities separated by commas and no spaces.
    * Example: ```MASSA_Algorithm -i <input_file>.sdf -o <output_file>.sdf -s pIC50,pMIC```.
    * Default = If not entered directly on the command line, it will be requested during algorithm execution.
* **Number of principal components in PCA**: ```-n``` or ```--number_of_PCs```.
    * Defines the number of principal components to reduce the dimensionality of variables related to biological, physicochemical and structural domains. If the value is a decimal between 0 and 1, the number of principal components is what justifies for (```<input number>```* 100)% of the variance. If the value is greater than 1, the number of PCs will be exactly the input integer, but PAY ATTENTION:

        1) If the number of PCs is an integer and equal to or greater than the number of physicochemical properties (7), the PCA step will be bypassed for this domain.
        2) The same for the biological domain.
        3) If the number of biological activities is less than 3, the PCA step will be bypassed for this domain.
    * Default = 0.85.
* **SVD solver parameter for PCA**: ```-v``` or ```--svd_solver_for_PCA```.
    * See the sklearn.decomposition.PCA topic on https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html for more info.
    * Default = full.
* **Extension of image files**: ```-t``` or ```--image_type```.
    * Extension of the image files that will be generated. Suggested = png or svg.
    * Default = png.
* **Font size for X-axis of dendrograms**: ```-d``` or ```--dendrogram_Xfont_size```.
    * Sets the font size on the x-axis of the dendrogram (molecule labels).
    * Default = 5.
* **Font size for X-axis of bar plots**: ```-x``` or ```--barplot_Xfont_size```.
    * Sets the font size on the x-axis of the bar plot (cluster labels).
    * Default = 12.
* **HCA linkage method**: ```-l``` or ```--linkage_method```.
    * The linkage criterion to use. The algorithm will merge the pairs of cluster that minimize this criterion.
    * Options = complete, single, ward, average, weighted, centroid, median. For more info, see the scipy.cluster.hierarchy.linkage topic on https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage.
    * Default = complete.
* **Enable Dendrogram plot**: ```-f``` or ```--dendrogram_plot```.
    * Defines whether or not dendrogram images will be generated.
	* Options = true (dendrogram will be generated), false (dendrogram will not be generated).
    * Default = true.

#### Command line help
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
