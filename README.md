# ISEtools

**ISEtools** is an R package for the statistical analysis of data from **ion-selective electrodes (ISEs)** using Bayesian methods.  
It is designed to make best-practice analysis of ISE calibration and experimental data accessible to non-specialists, while still allowing advanced users to customise models and workflows.

The package implements Bayesian models based on the **Nikolskii–Eisenman equation**, enabling robust estimation of sensor parameters, analyte activities, and limits of detection (LOD), including appropriate uncertainty quantification.

---

## Key Features

- **Bayesian analysis of ISE data** using OpenBUGS or JAGS
- Handles **non-linear sensor response**, including data near the limit of detection
- Supports:
  - Single ISEs or **arrays of multiple ISEs**
  - Calibration data only, or calibration + experimental samples
  - Experimental data collected via **basic** or **standard addition** methods
- Estimation of:
  - Model parameters (baseline, slope, interference parameter, noise)
  - **Limits of detection (LOD)** following IUPAC-consistent definitions
  - Analyte activities for experimental samples, with credible intervals
- **Substantial automation**: minimal Bayesian or BUGS expertise required
- Publication-ready summaries and plots via standard `print()`, `summary()`, and `plot()` methods

---

## Core Functions

ISEtools revolves around three main functions:

- `loadISEdata()`  
  Imports and validates calibration and (optional) experimental data from tab-delimited text files.

- `describeISE()`  
  Uses calibration data to estimate ISE model parameters and limits of detection.

- `analyseISE()`  
  Estimates unknown analyte activities in experimental samples, optionally combining information from multiple ISEs.

Each function returns a structured R object with associated `print`, `summary`, and `plot` methods.

---

## Installation

### From CRAN
```r
install.packages("ISEtools")
```

### From GitHub
### Development version (GitHub)
```r
# install.packages("remotes")
remotes::install_github("dillingh1/ISEtools")
```

---

## System Requirements

ISEtools relies on external Bayesian software:

- **R** (or RStudio)
- One of:
  - **OpenBUGS** (recommended; Windows/Linux)
  - **JAGS** (cross-platform; recommended for macOS)
- along with either of the corresponding R packages:
  - R2OpenBUGS
  - rjags   

---

## Documentation

- **Package vignette**:  
  ```r
  vignette("ISEtools")
  ```
- Function help pages:
  ```r
  ?loadISEdata
  ?describeISE
  ?analyseISE
  ```

The vignette contains detailed explanations of data formats, Bayesian models, examples, and advanced customisation options.

- **Published article describing use and primary citation for the software**:

> Dillingham, P.W., Alsaedi, B.S.O., Radu, A., & McGraw, C.M. (2019).  
> *Semi-Automated Data Analysis for Ion-Selective Electrodes and Arrays Using the R Package ISEtools.*  
> **Sensors**, 19(20), 4544.


---

## License

This package is released under the **GPL-3**.

---

## Contact

**Peter W. Dillingham**  
Department of Mathematics and Statistics  
University of Otago  
📧 peter.dillingham@otago.ac.nz

