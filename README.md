# How to create a multiModel OmicNavigator study

This repository contains an example of how to use the [OmicNavigator R
package][on-rpkg] to convert a publicly available multi-model omic analysis into 
a study package to be explored with the [OmicNavigator app][on-app]. For more 
details, please see the User's Guide attached to the [latest release][latest].

[on-rpkg]: https://github.com/abbvie-external/OmicNavigator
[on-app]: https://github.com/abbvie-external/OmicNavigatorWebApp
[latest]: https://github.com/abbvie-external/OmicNavigator/releases/latest

**Files:**

* [`setup.R`](./setup.R) - Installs the required R packages for the example

* [`data/`](./data/) - The input differential expression data (excel file) as 
downloaded from article ["The Global Phosphorylation Landscape of SARS-CoV-2 
Infection"][paper], published in CELL, Vol. 182 (3), pages 685-712.e19. 2020.

[paper]: https://www.sciencedirect.com/science/article/pii/S0092867420308114?via%3Dihub#app2

* [`build.R`](./build.R) - Builds the OmicNavigator study package from the
data found in [`data/`](./data/). Installs the study package and starts the web 
app.

## Run the code

Follow the steps below to install the dependencies and create the OmicNavigator 
study package. 

1. Install R package dependencies

    ```
    source("setup.R", local = new.env())
    ```

1. Create and install the OmicNavigator study package. This reads the analysis
results files in `data/`, converts them to an OmicNavigator study package,
installs the package, and starts the app.

    ```
    source("build.R")
    ```

## Acknowledgements

The example was adapted from the article ["The Global Phosphorylation Landscape 
of SARS-CoV-2 Infection"][paper], published in Cell, Vol. 182 (3), 
pages 685-712.e19. 2020.

If you use the data, please cite:

> Bouhaddou M, Memon D, Meyer B, et al.: The Global Phosphorylation Landscape 
> of SARS-CoV-2 Infection. Cell. 2020; 182(3): 685-712.e19.
https://www.sciencedirect.com/science/article/pii/S0092867420308114

