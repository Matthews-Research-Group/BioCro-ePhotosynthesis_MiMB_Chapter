# BioCro-ePhotosynthesis_MiMB_Chapter
This repository contains all scripts and data that are needed to repeat the examples shown in the chapter **Coupling Detailed Photosynthetic Kinetics to Crop Growth Models: The BioCro-ePhotosynthesis Framework** in the book **Plant Systems Biology**.

The folder structure is as follows,
- **models**: submodules including two models: BML-ePhotosynthesis (interface of the coupled BioCro-ePhotosynthesis) and ePhotosynthesis_C (the ePhotosynthesis C++ source).
- **parameterization**: scripts and data for calibrating and parameterizing the leaf-level ePhotosynthesis.
- **sanity_check**: scripts for logging and diagnosing the metabolites build-ups (convergence) when running the standalone ePhotosynthesis or the coupled BioCro-ePhotosynthesis model.
- **seasonal_simulation**: scripts for running and viewing seasonal simulations with the coupled model.
