[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13134608.svg)](https://doi.org/10.5281/zenodo.13134608)
[![DOI](https://img.shields.io/badge/Paper--DOI-10.5194%2Fegusphere--2024--2444-%23FFDE00.svg)](https://doi.org/10.5194/egusphere-2024-2444)  
[![output data regression test](https://github.com/opinner/Pinner_et_al_2024/actions/workflows/output_regression_tests.yml/badge.svg)](https://github.com/opinner/Pinner_et_al_2024/actions/workflows/output_regression_tests.yml)

# Analysis code to: [Pinner et al., 2024](https://doi.org/10.5194/egusphere-2024-2444)

This is the analysis code to the for open peer-review submitted publication  
[**Pinner et al., 2024:** *Internal-wave-induced dissipation rates in the Weddell Sea Bottom Water gravity current*](https://doi.org/10.5194/egusphere-2024-2444), 


## Disclaimer
- Although this code produces the results and figures to the accompanying paper, this repository still contains unused code snippets and unfinished documentation. 
The repository will be cleaned up further during the review process, while making sure that the output remains the same.
This may include the renaming of files or functions to better communicate their purpose.
- Comments or corrections to the code can be given on GitHub as issues or discussions or on EGUsphere as a community reviewer.  
- Note that figures created here can differ slightly from the submitted versions, as some post-processing (adjustements and labeling) were made with *Inkscape*. 

## Structure
The figures 1 to 6 are produced from data as follows:
```mermaid
flowchart TD
    bathymetry[(bathymetry)] ----> f01((f01: overview))
    CTD1[(CTD profiles)] --> id[[""eos80_legacy_gamma_n\nmatlab toolbox""]]
    id --> gamma["neutral density γⁿ"]
    gamma --> f01

    f02((f02: spectrum))
    TS[(velocity\ntime series)] ----> f02
    TS ----> f03((f03: flowfield))

    TS --> idemix[[wave energy/IDEMIX \nparameterization]]
    idemix --> epsidemix["ε_{IGW, IDEMIX}"]

    CTD2[(CTD profiles)]
    CTD2 --> fine[[finestructure]]
    fine --> epsfine["ε_{IGW, fine}"]
    f05(("f05: ε_{IGW}\ntransect"))
    epsfine --> f05
    epsidemix --> f05
    epsidemix --> f04(("f04: ε transect"))
    CTD2 --> Thorpe[[Thorpe scales]]
    Thorpe --> epstotal["ε_{total}"]
    epstotal --> f04

    f06(("f06: ε profile"))
    epstotal --> f06
    epsidemix --> f06
    epsfine  --> f06
```

* Data (1st level) is referenced from the `./data` directory (but not committed to repository)
* Methods (2nd level) are in the `./scripts` directory
* Methods results (3rd level) are in the respective `./scripts/<method_name>/method_results` directories
* Figures (4th level) are produced by the scripts in the `./results` directory
