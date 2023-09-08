# CIDMD_analysis
#### Authors: Jesi Lee

## General Description
* These scripts analyze CIDMD (collision-indueced dissociation molecular dynamics) simulations in TeraChem.



## Organization
For the scripts to work properly, 
* `Molecule/mol_info.in` is required (for details, see `CIDMD_setup/README.md` in the repository `github.com/jesilee/CIDMD_setup`.
*  The CIDMD trajectories should have been processed with `LearnReactions.py`, and this repository should be located at `Molecule/cid/CIDMD_analysis` and 
*  Each resulting subfolder `Molecule/cid/results/` will be created containing these resulting files:
  `cidmd.jdx`, `cidmd.msp`,
  `cidmd.PS_case1_RCID_vfac8-cut05norm.png`, `cidmd.PS_case1_RCID_vfac8-cut05.png`,
  `ar_vel_KE.png`,`sy_found.out`,
  `total_gputime.txt`, `pop_found.out`,
  `rxn_found.out`, `report_found.out`.



## Inputs
* `learn_rxn.log` files located in `Molecule/cid/calcs/XXX/gathered/` that resulted from post-processing CIDMD by running `LearnReactions.py` written by Prof. Lee-Ping Wang. His repository: 'https://github.com/leeping/nanoreactor' (not yet released)
* `mol_info.in` file located in `Molecule/`

  
## Outputs
* `results` folder under `Molecule/cid/` containing these files: `cidmd.jdx`, `cidmd.msp`, `cidmd.PS_case1_RCID_vfac8-cut05norm.png`, `cidmd.PS_case1_RCID_vfac8-cut05.png`, `ar_vel_KE.png`,`sy_found.out`, `total_gputime.txt`, `pop_found.out`, `rxn_found.out`, `report_found.out`.
  
* `cidmd.jdx`, `cidmd.msp`: CIDMD theoretical mass spectrum in JDX and MSP file format. Both file format are useful for the NIST programs or other mass spectra libraries.
* `cidmd.PS_case1_RCID_vfac8-cut05norm.png`, `cidmd.PS_case1_RCID_vfac8-cut05.png`: PNG format pictures of normalized and not normalized CIDMD theoretical mass spectra.
* `ar_vel_KE.png`: a PNG format picture of argon velocities as a function of collision time (usually less than 1 picosecond).
* `sy_found.out` : a survival yield value for molecular ion during CIDMD trajectories.
* `total_gputime.txt` : total gpu compute time for CIDMD in seconds and hours.
* `pop_found.out` : calculated charges of each fragments during CIDMD trajectories
* `rxn_found.out` : detected reactions by LearnReactions.py during CIDMD trajectories
* `report_found.out` : a summary of fragments and thier charges during CIDMD trajectories.


## Dependencies
* python 3.6
* numpy
* pandas
* matplotlib
* these files from this repo should be in the current working directory:
   `CIDMD_analysis.com`
   `CIDMD_analysis.py`
* `LearnReactions.py` written by Prof. Lee-Ping Wang. His repository: 'https://github.com/leeping/nanoreactor' (not yet released)
* `Molecule/cid/calcs/XXX/gathered/learn_rxn.log`
  (this would be the resulting output file from running `LearnReactions.py`)
* `Molecule/mol_info.in` (please see 'https://github.com/jesilee/CIDMD_setup/blob/main/README.md' for an example)

  
## Important notes
To start CIDMD, Folder organization is very important. Please see the instruction to set up CIDMD here:
'https://github.com/jesilee/CIDMD_setup/blob/main/README.md'

