# FLEXR
[alt text](https://github.com/TheFischerLab/FLEXR/tree/main/img/logo.png?raw=true)

`FLEXR` is a command line tool for building in alternative side chain conformations
determined from [Ringer](https://bl831.als.lbl.gov/ringer/ringer/Documentation/ringerManual.htm).

If you use this software, please cite:
[Placeholder]()

## Installation

`FLEXR` was tested on Intel and M1 Macs running macOS Monterey.
You will need the following tools:
1. git
2. [Phenix](https://phenix-online.org)
3. [Homebrew](https://brew.sh)

Once these are installed, you can:
1. Install dependencies (NumPy, SciPy, and Pandas):
```
pip install numpy
```
2. Clone the latest release of FLEXR:
```
git clone https://github.com/TheFischerLab/FLEXR.git
```
3. Install [Coot](https://pemsley.github.io/coot/blog/2022/06/05/coot-1-on-macos.html) (v1.0.5) (Emsley & Cowtan, 2004) via Hombrew using a formula developed by [Yoshitaka Moriwaki](https://github.com/YoshitakaMo):
```
brew install ./coot.rb --verbose --debug --keep-tmp
```
Helpful information for troubleshooting the Coot 1.0 installation can be found [here](https://github.com/pemsley/coot/issues/33).

4. Install Pandas to Coot packaged python.
Please note that the path on your computer might be slightly different:
```
/opt/homebrew/bin/python3.10 -m pip install pandas
```

## Usage examples

`FLEXR` relies on electron density measurements calculated by Ringer.
Ringer is packaged in the MMTBX module of the CCTBX library.
Maps suitable for Ringer can be produced with:
```
phenix.maps somepdb.pdb somesf.mtz
```
Ringer is then run with:
```
mmtbx.ringer somepdb.pdb some_map_coeffs.mtz
```
This produces a CSV file, `somepdb_ringer.csv`, which contains the raw electron
density measurements that `FLEXR` relies on.

### 1. Conformer detection

`FLEXR` contains two scripts, one for (1) conformer detection and one for (2) model building.
The first script is `flexr_find.py` which performs three main functions from
Ringer based electron density measurements. First, it performs peak detection, second
it assembles these peaks into possible rotamers, and third it tests these rotamers
against the ideal rotamer library (`rotamer_library_coot.csv`). It is run with:
```
python flexr_find.py -f somepdb_ringer.csv
```
Two options users will want to test are the electron density threshold `-t` for peak detection and
the geometry tolerance `-g` used for matching the ideal rotamer library:
```
python flexr_find.py -f somepdb_ringer.csv -t 0.35 -g 40
```
Parameter optimization can be performed using the `flexr_grid_optimization.sh` script in the `script` folder.
Rotamers slated for building are saved to `somepdb_ringer_alts.csv`.

Plots showing the electron density and detected peaks and be produced by setting the `-p` flag to True.

### 2. Model building

The second script is `flexr_build.py`, which takes the rotamers identified in the previous step
and builds them into a single conformer model using model building tools in Coot.
Please note that the path to Coot 1.0 might vary on your computer. This step is run with:
```
/opt/homebrew/Cellar/coot/1.0.05/bin/coot --script flexr_build.py file.somepdb
```

where `somepdb` is the name that matches both the input single conformer model file (.pdb) you want to build on and the prefix of the `_ringer_alts.csv` file containing the list of conformers that will be built.
For example, if your model name and `_ringer_alts.csv` prefix are 1ABC, use `file.1ABC`.

Another important setting is whether to add alternative conformations starting at the Cα atom or create an entirely new residue. The default is to create an entirely new residue, but this can be changed using:

```
/opt/homebrew/Cellar/coot/1.0.05/bin/coot --script flexr_build.py file.somepdb ca_or_all.0
```

### 3. Refinement

Output multi-conformer models are named `somepdb_newconfs.pdb` and can be refined using the program of the user’s choice.
The Phenix PDB editing tool, `phenix.pdbtools`, can be used to filter out rotamers with certain occupancies.
For example, to remove rotamers with occupancy=0 use:
```
phenix.pdbtools somepdb_newconfs_refined.pdb remove="occupancy=0”
```
