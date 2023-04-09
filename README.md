### DiseasePredictor

This package was developed as part of the exam of the Scientific Programming course, taken in the 
"Bioinformatics for Computational Genomics" M.Sc of Politecnico di Milano.

#### Introduction
This package implements functionalities for quantifying the likelihood of some genes
(the <i>candidate genes</i>) of being related to a particular disease. This computation is 
based on a protein-protein interaction network, coming from the STRING database (https://string-db.org/),
and representing gene products' molecular interactions, and a list of <i>seed genes</i>. Each candidate gene
is assigned a score based on the average distance from each seed gene node on the network.
The distance metric I decided to adopt is the <i>random walk with restart</i>, whose formulation in this context
is well described in the study <i>"Köhler, Sebastian et al. 'Walking the interactome for prioritization 
of candidate disease genes.' American journal of human genetics vol. 82,4 (2008)"</i>.

#### Project structure
The code is structured as follows:
<ul>
    <li>The source package, named <i>disease_predictor</i>, and
    it contains all the classes and functions needed for building a
    representation of the protein-protein interaction network, computing the
    scores, ranking genes and plotting the scores.
        <ul>
        <li>The <b>DiseasePredictor</b> class is the engine of the package.
            When initialized with a specific configuration, it will build
            the internal representation of the PPI network, and expose methods 
            for retrieving the network, ranking candidate genes based on it, and
            obtain a simple plot of the scores.<br/>
            The configuration is given by instantiating the <b>DiseasePredictorConfig</b>
            class, that has to contain the network's filename and a protein-gene mapping's filename.
            They must be tabular files with a heading, thus some characters must also be provided
            as separators for the parsing.<br/>
            The network file must follow the STRING format, and thus the first three
            columns must contain the two interacting proteins and the score for the interaction.<br/>
            The gene-protein mapping will simply need to contain, for each protein ID in the network, the 
            corresponding official gene symbol.<br/>
            Other optional parameters to pass to the configuration are a threshold for the interaction scores,
            and the number of steps and the restart probability for the random walk implementation.
            </li>
        </ul>
    </li>
    <li>The <i>tests</i> package contains some unit tests written
    for the main functions, with some simple test cases designed for
    ensure the methods return the expected results.
    </li>
</ul>

To run the tests, execute
```shell
python3 -m pytest test/
```

#### How to use the package
You can find the already built package, ready to 
be installed and used, inside the <i>demo</i> directory.
To re-build the package, it is sufficient to run
```shell
python3 -m build
```
in the source code root directory. The package will be
generated inside the <i>dist</i> directory.

Once you have the .tar.gz package, you can simply install it by
running
```shell
pip install /path/to/the/pckg/disease_predictor-1.0.0.tar.gz
```

And then import it inside your python source code with
```python
import disease_predictor
```

#### Demo
Along with the built package, in the <i>demo</i> directory you can also find a simple 
demo script (<i>demo_disease_predictor.py</i>) I wrote to
provide a usecase of the package.
I decided to test the package with genes associated with two diseases,
diabetes and depression. I downloaded the gene lists after performing two phenotypes searches on
https://www.ncbi.nlm.nih.gov/gap/PheGenI/, one for each of the two pathologies. 
I built the network based on a STRING db PPI network
and set as seed genes 100 diabetes genes, then assessed whether diabetes candidate genes
showcased higher scores with respect to depression candidate genes.
The plot, that can be also found inside the directory (<i>plot_1673364971.png</i>),
shows how many diabetes genes seem to indeed have higher scores.
