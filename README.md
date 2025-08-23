[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](
https://colab.research.google.com/github/fowler-lab/tuberculosis-RNAP-subpopulations/blob/main/analysis_refactored.ipynb
)

# tuberculosis-RNAP-subpopulations
This repository contains a juypter-notebook, attendant code and data allowing all analysis and (nearly all) figures in the below preprint to be reproduced:

Brunner VM, Fowler PW
**Subpopulations in clinical samples of M. tuberculosis can give rise to rifampicin resistance and shed light on how resistance is acquired**
bioRxiv preprint [doi:10.1101/2025.04.09.647945](https://doi.org/10.1101/2025.04.09.647945)

## Instructions

### To run in browser
The easiest way is to click the blue `Open in Colab` button at the top left of this README when viewed on github.com. This will open the notebook in your browser window and by clicking the "play" buttons in the grey code windows you can step through the code, including drawing the graphs etc.

### To run locally
You need have installed a number of standard Python packages -- many of these are common (e.g. `numpy`) so you may already have a version installed. To install these locally using `pip` issue

```
$ git clone git@github.com:fowler-lab/tuberculosis-RNAP-subpopulations.git
$ cd tuberculosis-RNAP-subpopulations/
$ pip install -r requirements.txt
```

Then you will need to open the notebook using your IDE of choice. For example, to open in VS Code

```
$ code -n notebook.ipynb
```

And you can step through the different code cells as above using the "play" buttons.

Running the entire notebook should take no more than a few minutes.

## Source of data
As mentioned in the manuscript, the data tables included in this repository are extracts from [v3.0.0](https://doi.org/10.5281/zenodo.16041005) of the CRyPTIC Consortium Dataset.

