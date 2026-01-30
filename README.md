# MDS-project
pipeline for plant TIR domains (ham8)
This repository/workspace contains a reproducible pipeline to:
1) run FrustratometeR configurational frustration on a large set of TIR-domain PDBs,  
2) aggregate contact-level outputs into per-residue metrics,  
3) generate cross-model plots (by absolute residue number and by domain-relative `tir_pos`),  
4) export residue attributes for UCSF ChimeraX (`.defattr`) for colouring and figure export.

## 1. HPC environment (ham8)

### 1.1 Load modules (interactive)
```bash
module purge
module load r/4.4.1
module load python/3.9.9    # or python/3.10.8 if preferred
module load chimerax/1.7.1  # optional, only for visualisation
1.2 Personal R library path 
Set a per-user R library to avoid permissions issues:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
mkdir -p "$R_LIBS_USER"
1.3 Python venv 
python -m venv "$HOME/venvs/tir_frustra_py"
source "$HOME/venvs/tir_frustra_py/bin/activate"
pip install -U pip pandas numpy

2. Directory structure
Expected structure (one folder per species):
~/tir_frustra/
├── species/
│   ├── arabidopsis/
│   ├── Nicotiana_tabacum/
│   ├── Solanum_lycopersicum/
│   ├── Vitis_vinifera/
│   │   ├── pdbs/
│   │   ├── pdb_list.txt
│   │   ├── results/
│   │   │   └── <MODEL_ID>/
│   │   │       ├── <MODEL_ID>.rds
│   │   │       ├── DONE
│   │   │       └── <MODEL_ID>.done/
│   │   │           └── FrustrationData/
│   │   │               └── <MODEL_ID>.pdb_configurational
│   │   └── logs/
└── postprocess/
Model naming convention
Input model IDs encode domain boundaries:
AF-<accession>-TIR_<start>-<end>_buffer10
Example:
AF-A0A438BLY1-TIR_18-158_buffer10
3. Installing FrustratometeR (ham8)
In R (ensure R_LIBS_USER is set as above):
install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("protein-tools/FrustratometeR")  # repo name may vary in your setup
installation may complain about system libraries (e.g., magick):
•	install the dependency if supported on ham8, or
•	configure FrustratometeR to run without optional plotting/image steps (recommended for HPC batch).
Sanity check:
library(FrustratometeR)
packageVersion("FrustratometeR")
4. Running frustration test
From a species folder (example: arabidopsis):
cd ~/tir_frustra/species/arabidopsis
PDB=$(sed -n '1p' pdb_list.txt)
OUTDIR="$PWD/results/TEST_SINGLE"

Rscript "$PWD/run_one_frustra.R" "$PDB" "$OUTDIR"
ls -la "$OUTDIR"
Expected outputs include:
•	<MODEL_ID>.rds
•	<MODEL_ID>.done/FrustrationData/<MODEL_ID>.pdb_configurational
•	DONE marker file
5. Batch processing on SLURM
5.1 Prepare pdb_list.txt
pdb_list.txt should contain one model ID per line (matching files in pdbs/).
Example:
AF-A0A438BLY1-TIR_18-158_buffer10
AF-A0A3Q7JE37-TIR_10-177_buffer10
...
5.2 Example SLURM array script
Save as run_frustra_array.slurm inside a species folder:
#!/bin/bash
#SBATCH --job-name=frustra
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-1000
#SBATCH --output=logs/frustra_%A_%a.out
#SBATCH --error=logs/frustra_%A_%a.err

module purge
module load r/4.4.1
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"

SPECIES_DIR="$(pwd)"
MODEL_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pdb_list.txt)

# skip empty lines
if [ -z "$MODEL_ID" ]; then
  echo "No MODEL_ID for index ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

PDB_PATH="${SPECIES_DIR}/pdbs/${MODEL_ID}.pdb"
OUTDIR="${SPECIES_DIR}/results/${MODEL_ID}"

# idempotency: skip if DONE exists
if [ -f "${OUTDIR}/DONE" ]; then
  echo "Already completed: ${MODEL_ID}"
  exit 0
fi

mkdir -p "${OUTDIR}"
Rscript "${SPECIES_DIR}/run_one_frustra.R" "${PDB_PATH}" "${OUTDIR}"
Submit:
mkdir -p logs
sbatch run_frustra_array.slurm
Post-processing (per-residue aggregation + plots)
6.1 Inputs used as “authoritative”
We treat the contact-level output as authoritative:
results/<MODEL_ID>/<MODEL_ID>.done/FrustrationData/<MODEL_ID>.pdb_configurational
6.2 Per-residue metrics
Each contact contributes to both residues. For each residue we compute:
•	frst_index (mean frustration index over contacts)
•	prop_highly (proportion of contacts labelled highly frustrated)
•	prop_neutral
•	prop_minimally
•	n_contacts
Per species output:
•	<species>_per_residue_from_configurational.csv
All species combined:
•	all_species_per_residue_from_configurational.csv
6.3 Domain-relative mapping (tir_pos)
To compare across different domain boundaries, compute:
tir_pos = resno - tir_start + 1
Outputs:
•	freq_high_by_resno.csv (and png)
•	freq_high_by_tirpos.csv (and png)
•	mean_frst_index_by_tirpos.csv (and png)
BB-loop definition 
BB-loop is defined structurally as the loop connecting βB to αB.
Conservative window used across all models:
•	tir_pos 34–49
8. ChimeraX visualisation
8.1 Important format note
ChimeraX does not read .attr files (UCSF Chimera legacy). Use .defattr.
Minimal header example:
attribute: prop_highly
recipient: residues
Data lines must be:
<TAB>resno<TAB>value
Example line:
    58    0.125
8.2 Load attributes
In ChimeraX:
open AF-A0A438BLY1-TIR_18-158_buffer10.pdb
defattr prop_highly.defattr
8.3 Colour mapping (low → high)
Palette left→right maps low→high:
color white
color byattribute prop_highly palette white:red range 0,1
8.4 Highlight BB-loop (example selection)
Replace residue numbers with the model’s absolute resno range:
select :51-66
color yellow cartoons sel
label :58 text "BB-loop"
9. Reproducibility notes
•	All analyses are apo monomer baseline unless explicitly extended.
•	Frustration is reported as relative patterns (enrichment), not absolute energetics.
•	Use DONE markers to ensure reruns are idempotent.
•	Keep module versions recorded:
o	r/4.4.1
o	chimerax/1.7.1
o	python/3.9.9 (optional)
10. Troubleshooting
“FrustratometeR installation fails (system deps)”
•	Ensure R_LIBS_USER is set
•	Try reinstalling in a clean session
•	If failures mention optional packages (e.g., image/plotting), disable optional steps in your run scripts for HPC.
“ChimeraX didn’t colour residues”
•	Confirm you used .defattr, not .attr
•	Confirm the data lines are TAB-separated
•	Confirm residue numbering matches the PDB resno values
11. Quick start (copy/paste)
module purge
module load r/4.4.1
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"

cd ~/tir_frustra/species/arabidopsis
mkdir -p logs
sbatch run_frustra_array
