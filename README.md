
This repository contains the code used in a Master's thesis to initialize a forest model (iLand) using LiDAR data.

## Structure

The workflow is designed to preprocess input files for each AOI (Area of Interest) separately and store them safely. These files can then be copied into an iLand project to run simulations. After the simulation, the output can be copied back to the AOI folder.

### Folders

* **allometries**: Contains code to test different allometric models and extract species-specific coefficients.
* **init**: Contains everything needed for initialization, including creating input files, an optional spin-up routine, and updating project files.
* **setups**: Contains project files for two different setups and climate scenarios.

### Scripts

* **copy_model_input.sh**: Moves initialization files to the selected iLand project.
* **copy_model_output.sh**: Moves output files from the iLand project to the AOI directory.

## Workflow

1. **Allometries**

   * Use `allometric_models.sh` in the `allometries` folder to test different regression models and extract species-specific coefficients for predicting DBH from tree height.

2. **Preprocessing**

   * Prepare all required input files for an AOI using the scripts in the `init` folder.
   * Run `run_prep.sh` to generate the main input files (calls `init_iland.R` and `climate_data.R`).

3. **Project initialization**

   * Copy all prepared input files into the selected iLand project directory and update project file using `copy_model_input.sh` (calls `copy_input.sh` and `parameter_xml.R`).

4. **Simulation**

   * Run iLand using the project files from the `setups` folder.

5. **Collect output**

   * Use `copy_model_output.sh` to copy the simulation results back into the AOI folder.
   * Store and organize outputs for later analysis.

6. **Optional: Spin-up routine**

   * Run an iLand simulation with regeneration enabled using a historic climate setup for 500 years (e.g., `setups/withSpinup/projectFile_historic.xml`).
   * Use the resulting long-term output as input for the spin-up script: `init/run_spinup.sh` (calls `spinup_treeinit.R`).
   * The spin-up script generates a new tree inventory.
   * Replace the initial tree inventory with the spin-up inventory and repeat steps 3–5 for final simulations.
