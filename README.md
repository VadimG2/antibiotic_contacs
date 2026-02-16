# Molecular Dynamics Simulation of Peptide-Antibiotic Interactions (GROMACS)

This repository contains a set of scripts for automated molecular dynamics (MD) simulation of peptide–antibiotic complexes using GROMACS, followed by analysis of the results.

## Pipeline Overview

For every peptide–ligand pair found in the input directories, the pipeline performs the following steps:

### Peptide Preparation
- Generate GROMACS topology (`.top`), coordinates (`.gro`), and position restraint file (`posre_*.itp`) from the peptide PDB file using `gmx pdb2gmx`.
- **Manual check/edit** of the `.top` file is required to ensure cyclization of the peptide (if applicable).
- Add the directive `#ifdef POSRES` to the `.top` file to enable the use of the position restraint file.

### Peptide Optimization (Annealing)
- Run a short vacuum MD simulation (annealing) to relax the initial peptide structure.
- Extract the last frame as the optimised structure (`_optimized.gro`).

### Ligand Parameterisation
- Convert the ligand PDB file to MOL2 format using Open Babel.
- Assign GAFF2 atom types and compute BCC charges with `antechamber` (from AmberTools).
- Find missing parameters with `parmchk2`.
- Build AMBER topology (`.prmtop`) and coordinates (`.inpcrd`) using `tleap`.
- Convert AMBER files to GROMACS format (`.itp`, `_GMX.gro`) with ACPYPE.
- Extract GAFF2 atom type definitions (`atomtypes_*.itp`) and position restraint file (`posre_*.itp`) from ACPYPE intermediate files.
- Correct the molecule name and the path to the position restraint file inside the ligand `.itp` file.

### System Assembly
- Create a unique simulation directory (`output/PEPTIDE+LIGAND_Sim`).
- Copy the necessary peptide and ligand files into the simulation directory.
- Merge peptide and ligand coordinates (`insert-molecules`).
- Create the simulation box (`editconf`).
- Automatically generate the `system.top` file:
  - Include required force field files (`amber14sb.ff`), GAFF2 atom types (`gaff2_atomtypes.itp`), ligand (`ligand.itp`), water (`tip3p.itp`), and ions (`ions.itp`).
  - Include the peptide `moleculetype` definition (with its `PosRe`).
  - Add the `[ system ]` and `[ molecules ]` sections (with peptide and ligand).
- Solvate the system with water (`solvate`).
- Add ions to neutralise the system and reach the desired salt concentration (`grompp`, `genion`).

### MD Simulation (executed by a separate script)
- Move the assembled system to the `simulations/` directory.
- Create an index file (`make_ndx`) with groups `Protein_Ligand` and `SOL_Ion`.
- Copy and adapt MDP templates (for minim, nvt, npt, prod) with the correct temperature/pressure coupling groups.
- Run the stages sequentially:
  - Energy minimisation (`grompp`, `mdrun`).
  - NVT equilibration (`grompp`, `mdrun` with position restraints).
  - NPT equilibration (`grompp`, `mdrun` with position restraints).
  - Production MD (`grompp`, `mdrun` without position restraints).

### Analysis (executed by a separate script)
- Generate contact maps between peptide residues and the ligand from the production trajectory using Python (MDAnalysis, matplotlib).

## Requirements

- **Conda**: Miniconda or Anaconda installed.
- **System dependencies**:
  - NVIDIA CUDA Toolkit (version compatible with GROMACS, e.g., 11.8 or newer)
  - Appropriate NVIDIA drivers
  - `wget` and `tar` utilities (usually present on Linux/WSL)
  - `bc` utility (usually present)
- **Git**: For cloning the repository.

## Installation and Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/VadimG2/antibiotic_contacts.git
   cd antibiotic_contacts
   ```

2. **Create and activate the Conda environment**:
   ```bash
   conda env create -f env.yaml
   conda activate gromacs_md_env   # use the environment name specified in env.yaml
   ```
   The `env.yaml` file lists all required Conda packages (GROMACS with CUDA, AmberTools, ACPYPE, Open Babel, Python libraries, etc.).

3. **Install the AMBER14SB force field**:
   The GROMACS package from Conda Forge may not include all force fields. This step downloads and installs AMBER14SB into the created environment.
   ```bash
   wget -O ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz https://ftp.gromacs.org/contrib/forcefields/amber14sb.ff.tar.gz
   tar -xvf ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz -C ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/
   rm ~/miniconda3/envs/gromacs_md_env/share/gromacs/top/amber14sb.ff.tar.gz
   ```
   *(Adjust the path if your Conda installation is not in `~/miniconda3`)*

4. **Prepare input files**:
   - Place your peptide PDB files in `input/protein/`.
   - Place your ligand (antibiotic) PDB files in `input/ligand/`.

5. **Clean previous results** (if you need to restart from scratch):
   ```bash
   rm -rf simulations/* output/*
   # This command deletes all previous results – use with caution!
   ```

6. **Run the main pipeline script**:
   ```bash
   bash query.sh   # or whatever your main script is named
   ```

## Directory Structure

```
.
├── input/
│   ├── protein/          # Raw peptide PDB files
│   └── ligand/           # Raw ligand PDB files
├── output/
│   ├── PEPTIDE_NAME/          # Prepared peptide files (.gro, .top, _optimized.gro, posre_*.itp)
│   ├── LIGAND_NAME/           # Prepared ligand files (.itp, _GMX.gro, atomtypes_*.itp, posre_*.itp)
│   └── PEPTIDE+LIGAND_Sim/    # Assembled system (before moving to simulations/)
├── simulations/
│   └── PEPTIDE+LIGAND/        # Working directory for a specific simulation (.tpr, .log, .xtc, .edr, index.ndx, …)
├── scripts/                    # (Recommended) Bash pipeline scripts (prepare_peptides.sh, prepare_ligands.sh, assemble_systems.sh, run_simulations.sh, analyze_contacts.py)
├── mdp_files/                  # (Recommended) Template .mdp files
├── env.yaml                    # Conda environment specification
└── README.md                   # This file
```

## Important Notes

- **File naming convention**: All peptide and ligand files should follow the format `peptX` and `ligX` (e.g., `pept1`, `lig1`), but they do not have to be numbers – patterns like `pept_something` and `lig_something` are also acceptable.
- The main pipeline script is `query.sh`. After successful system assembly, you will need to run the simulation script separately (see below).
- The pipeline assumes a **neutral or charged system**; ion concentration can be adjusted in the ion addition step.

## Running Simulations (Step 6)

After the system assembly script (e.g., `step5_...`) completes successfully and creates folders named `*_Sim` (which are then moved into `simulations/`), you need to launch the simulation script (e.g., `run_simulation_pipeline.sh v16` as discussed).

Before running `run_simulation_pipeline.sh`, ensure:

- The assembled system folders are inside the `simulations/` directory.
- Template MDP files (e.g., `step6.*.mdp`) are located in the project root directory.
- You are inside the activated Conda environment (`gromacs_md_env`).

Then execute:
```bash
bash run_simulation_pipeline.sh   # Runs minimisation, equilibration, and production
```

