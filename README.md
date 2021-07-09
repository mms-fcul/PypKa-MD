# PypKa-MD

PypKa + MD = constant-pH molecular dynamics

Implementation of the stochastic titration method <sup>1</sup>

[1] Baptista *et al.*, J. Chem. Phys. 117, 4184 (2002) DOI: 10.1063/1.1497164

## Installation

```
git clone https://github.com/mms-fcul/cphmd.git .
```

## Dependencies

Both PypKa and GROMACS are required to be installed in the system.

- PypKa >= 2.3.0
- GROMACS >=5.1.5

## Usage

PypKa-MD can be used by executing the pypkamd folder with a python3, giving as an argument a modified GROMACS .mdp input file to include Constant-pH specific variables.

```
python3 pypkamd System.mdp
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