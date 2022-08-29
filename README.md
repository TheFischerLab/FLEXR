# Ringer-refine
Multi-conformer modeling of crystallographic protein structures from electron density

Ringer-refine is a command line tool for building in alternative side chain conformations
determined from [Ringer]](https://bl831.als.lbl.gov/ringer/ringer/Documentation/ringerManual.htm).

If you use this software, please cite:
[Placeholder]()

## Installation

You will need the following tools:
1. git
2. [Phenix](https://phenix-online.org)
3. [Homebrew](https://brew.sh)

Once these are installed, you can:
1. Install dependencies (NumPy, SciPy, and Pandas):
```
pip install numpy
```
2. Clone the latest release of Ringer-refine:
```
git clone https://github.com/TheFischerLab/ringer-refine.git
```
3. Install Coot (v1.0.) (Emsley & Cowtan, 2004) via Hombrew using a formula developed by Yoshitaka Moriwaki:
```
wget https://raw.githubusercontent.com/YoshitakaMo/homebrew-bio/coot/Formula/coot.rb -O coot.rb
brew install ./coot.rb
```
4. Install Pandas to Coot packaged python:
```
/opt/homebrew/bin/python3.10 -m pip install pandas
```

## Usage examples

Ringer-refine relies on electron density measurements calculated by Ringer.
Ringer is packaged in the MMTBX module of the CCBTX library.
Maps suitable for Ringer can be produced with:
`phenix.maps somepdb.pdb somesf.mtz`
Ringer is then run with:
`mmtbx.ringer somepdb.pdb some_map_coeffs.mtz`
This produces a CSV file, `somepdb_ringer.csv`, which contains the raw electron
density measurements that Ringer-refine relies on.
### 1. Conformer detection
`Ringer-refine` contains several scripts, one for each step of model building.
The first script is `ringer_refine.py` which performs three main functions from
Ringer based electron density measurements. First, it performs peak detection, second
it assembles these peaks into possible rotamers and third it tests these rotamers
against the ideal rotamer library. It is run with:
`python ringer_refine.py -f somepdb_ringer.csv`
Two options users will want to test are the electron density threshold `-t` for peak detection and
the geometry tolerance `-g` used for matching the ideal rotamer library:
`python ringer_refine.py -f somepdb_ringer.csv -t 0.35 -g 40`
Rotamers slated for building are saved to somepdb_ringer_alts.csv. Plots showing the electron density and detected peaks and be produced by setting the `-p` flag to True.
### 2. Model building
The second script is `coot_ringer_build.py`, which takes the rotamers identified in the previous step
and builds them into a single conformer model using model building tools in Coot with:
`/opt/homebrew/Cellar/coot/1.0.05/bin/coot --script coot_ringer_build.py file.somepdb`
where ‘somepdb’ is name that matches both the input single conformer model file (.pdb) you want to build on and the prefix of the _ringer_alts.csv file containing the list of conformers that will be built.
Another important setting is whether to add alternative conformations starting at the Cα atom or create an entirely new residue. The default is to create an entirely new residue, but this can be changed using:
`/opt/homebrew/Cellar/coot/1.0.05/bin/coot --script coot_ringer_build.py file.somepdb ca_or_all.0`
Output models are name `somepdb_newconfs.pdb` and can be refined using the program of the user’s choice.
