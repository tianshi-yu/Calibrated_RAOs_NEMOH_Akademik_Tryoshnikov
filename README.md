# Paper Supplement: Calibrated RAOs for Akademik Tryoshnikov

## Introduction
Supplementary materials for the paper:

[T. Yu, Z. Wang, F. Nelli, Y. Tan, G. Ducrozet, and A. Toffoli, “Assessing the ship motion prediction capabilities of the open-source model NEMOH against field observations,” Ocean Eng., vol. 344, Art. no. 123580, Jan. 2026, doi: 10.1016/j.oceaneng.2025.123580.](https://www.sciencedirect.com/science/article/pii/S0029801825032627) 
 
This repository contains response amplitude operators (RAOs) and hydrodynamic parameters for the research vessel Akademik Tryoshnikov, computed using the open-source software NEMOH and calibrated against field observations during Antarctic Circumnavigation Expedition (ACE). For details please refer to the paper.


## Files
```text
├── main.m                      # Example: load and recompute calibrated RAOs
├── computeRAO.m                # Compute RAOs using calibrated parameters
├── calIkedaRollDamping.m       # Calculate roll damping using Ikeda's method
├── CalibratedRAO.mat           # Calibrated complex RAO arrays
└── CalibratedParameters.mat    # Calibrated/original parameters and speed/frequency/direction grids
```


## Usage
* Reuquirement: MATLAB

* Please refer to [main.m](main.m) for usage.

* We would appreciate your acknowledgement by citing our paper:
   
    ```    
   @article{YU2026123580,
    title = {Assessing the ship motion prediction capabilities of the open-source model NEMOH against field observations},
    journal = {Ocean Engineering},
    volume = {344},
    pages = {123580},
    year = {2026},
    issn = {0029-8018},
    doi = {https://doi.org/10.1016/j.oceaneng.2025.123580},
    url = {https://www.sciencedirect.com/science/article/pii/S0029801825032627},
    author = {Tianshi Yu and Ziyue Wang and Filippo Nelli and Ying Tan and Guillaume Ducrozet and Alessandro Toffoli}
    }
    ```
## Related Publication
NEMOH:
[LHEEA, 2022. NEMOH software v3.0, ecole centrale nantes, France. https://gitlab.com/lheea/Nemoh](https://gitlab.com/lheea/Nemoh)

Vessel motion data:
[Alberello, A., Bennetts, L., Toffoli, A., 2020a. [Database] Antarctic circumnavigation expedition 2017: motion sensor and GPS data . https://doi.org/10.4225/15/5A178EF0E5156](https://doi.org/10.4225/15/5A178EF0E5156)

WAMOS data:
[Alberello, A., Bennetts, L., Toffoli, A., 2020b. [Database] Antarctic circumnavigation expedition 2017: WAMOS data product . https://doi.org/10.26179/5E9D038C396F2](https://doi.org/10.26179/5E9D038C396F2)
## License

[![License: CC BY 4.0](https://licensebuttons.net/l/by/4.0/80x15.png)](https://creativecommons.org/licenses/by/4.0/) CC BY 4.0 

## Contact
Dr Tianshi Yu, tians.yu7@gmail.com
