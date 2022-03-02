# PypKa-MD

PypKa + MD = constant-pH molecular dynamics

Implementation of the stochastic titration method <sup>1</sup>

[1] Baptista *et al.*, J. Chem. Phys. 117, 4184 (2002) DOI: 10.1063/1.1497164

## Installation

```
python3 -m pip install pypkamd
```

## Dependencies

Both PypKa and GROMACS are required to be installed in the system.

- PypKa >= 2.7.1
- GROMACS >=5.1.5

When running in pKAI-MD mode there are extra dependencies:

- pege >= 1.1.1
- torch_geometric >= 2.0.0

Please refer to the installation guide of [torch geometric](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html) to install the proper version in accordance to your CUDA and OS.

```
python3 -m pip install pege
# EXAMPLE FOR CUDA10.2
# python3 -m pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.10.0+cu102.html
```


## Usage

Upon installation a PypKa-MD executable should have been added to your bin. You may call it directly giving as an argument a modified GROMACS .mdp input file to include Constant-pH specific variables.

```
pypkamd System.mdp
```

In case the executable as not been added to your bin, you may use:

```
python3 -m pypkamd System.mdp
```

You may find an example .mdp file in /utils/cphmd.mdp. 

```
; GROin = system_000.gro     ; input GRO file
; TOPin = system_000.top     ; input TOP file
; DATin = fixgro.dat         ; input DAT file (to be removed)
; NDXin = system.ndx         ; input NDX file
; sysname = system_001       ; output files root name
; sites = all                ; to be titrating sites in the protein
; titrating_group = Protein  ; index group of the protein
; nCycles = 50               ; number of CpHMD cycles
                            ;; total simulation time = nCycles * tau_prot
                            ;; 1ns = 50 * 20ps
; nCPUs = 4                  ; number of CPUs to be used
; pH = 7.0                   ; pH value of the protonation states sampling
; ionicstr = 0.1             ; ionic strength used in PB
; GroDIR="/gromacs/gromacs-5.1.5_pH_I/bin/" ; GROMACS bin path
```