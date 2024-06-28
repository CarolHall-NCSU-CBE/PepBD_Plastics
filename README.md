# Description
PepBD_Plastics is a Fortran program that discovers polypeptides predicted to have affinity for common plastics using Monte Carlo sampling and atomistic biophysical modeling. 

# Required input files
All should be placed in the directory for running PepBD, unless otherwise noted. Examples are provided in the repository
- compiled executable
- input.txt (see below for details)
- sys.pdb (see below for details)
- the lib folder, placed at the location ~/PepBD/
- A directory named `pdbfiles` 

# Installation and Compilation
Requirements: Intel ifort compiler (tested on versions 2017 - 2024)

  1. Navigate to the directory with PepBD_Plastics.f90
  2. Determine the total number of residues, `num_residues` in sys.pdb
  3. Replace `gnum=0` in with `gnum=num_residues` in  PepBD_Plastics.f90 file.
  4. Run the command `ifort -O2 -o main PepBD_Plastics.f90`

# Running PepBD_Plastics

Run from a terminal using `./main`. Note that PepBD_Plastics is a serial program.

The runtime depends on the system size and number of design steps, but normally takes a day to a week to complete.

# Output of PepBD_Plastics
Examples of each output are provided in the repository
   1. energyprofile.txt: score and energy components at each PepBD step
   2. energydetails.txt: sequence, energy, and move type at each PepBD step
   3. output.txt: Details of the PepBD run
   4. Generated pdbs: stored in the `pdbfiles` directory
   5. rmsd.txt: the rmsd relative to the starting state at each step of design. This is disabled in PepBD_Plastics.f90, so 0 is output at every step.
 
# Details on required files
# input.txt
The file must be named `input.txt`. The following variables that can be specified 
- PDBFILE (required): the name of the pdb file containing the peptide-receptor complex
- RESTARTPDB (required): the name of the pdb file if restarting PEPBD halfway through design. If starting at beginning of design, then name should match PDBFILE
- RECALCUSWITCH (required): Specifies if design is starting at beginning or being restarted from previous run. Set to 0 if starting from beginning, 1 otherwise. If set to 1, then the coordinates in "RESTARTPDB" will be used as the starting point for design
- STEP_START (required): Beginning step number of PepBD design. Set to 0 if starting design from beginning. Not used if SA_FLAG is set to 1, but still required.
- STEP_END (required): Final step number of PepBD design. Not used if SA_FLAG is set to 1, but still required.
- SITENUM (required): Number of residues in the peptide.
- RECEPTOR_NAME (required): Type of receptor. Options are nucleic, peptide, or other
- BACKBONE_CHANGE_PROB (required): Probability of making a change to peptide backbone for a given PepBD step.
- SEQ_CHANGE_PROB (required): Probability of making a change to the peptide sequence for a given PepBD step.
- SCMF_SWITCH (default=0.8): For a peptide sequence change, the probability that the change is changing an amino acid rather than switching two amino acids
- KT_SEQUENCE (required): The reference temperature for peptide sequence changes used in the Metropolis acceptance criterion.
- KT_BACKBONE_RATIO (required): Ratio of reference temperature for backbone changes to reference temperature fo sequence changes
- pH (required): pH of system. Hydration state of peptide and receptor should match the system pH
- RANSEED (required): seed to random number generator
- MIN_RMSD (required): Minimum RMSD limit for peptide backbone
- MAX_RMSD (required): Maximum RMSD limit for peptide backbone
- NUM_GLY_MIN (required): Lower limit on number of glycine that must be in peptide sequence
- NUM_GLY_MAX (required): Upper limit on number of glycine that must be in peptide sequence
- NUM_PHOBIC_MIN (required): Lower limit on number of hydrophobic residues that must be in peptide sequence
- NUM_PHOBIC_MAX (required): Upper limit on number of hydrophobic residues that must be in peptide sequence
- NUM_PHIL_MIN (required): Lower limit on number of hydrophilic residues that must be in peptide sequence
- NUM_PHIL_MAX (required): Upper limit on number of hydrophilic residues that must be in peptide sequence
- NUM_NEG_MIN (required): Lower limit on number of anionic residues that must be in peptide sequence
- NUM_NEG_MAX (required): Upper limit on number of anionic residues that must be in peptide sequence
- NUM_POS_MIN (required): Lower limit on number of cationic residues that must be in peptide sequence
- NUM_POS_MAX (required): Upper limit on number of cationic residues that must be in peptide sequence
- NUM_OTHER_MIN (required): Lower limit on number of other residues that must be in peptide sequence
- NUM_OTHER_MAX (required): Upper limit on number of other residues that must be in peptide sequence
- TRP_LIMIT (default=3): Upper limit on number of tryptophan in the peptide sequence
- NUM_HELIX (required): Number of helices in the peptide
- STEP_BEFORE_CONF_CHANGE (default=50): Number of PepBD steps to take before changes to the peptide backbone are allowed
- PDB_GEN_FREQ  (default=50): Number of PepBD steps taken before a pdb file of system is changed. Note that pdb files are also saved whenever the PepBD score is the lowest out of all sequences/backbones sampled
- PEP_WEIGHT_FACTOR (default=0.01): Factor that multiplies peptide-peptide interactions when calculating the PepBD score
- INITIAL_RANDOMIZE (default=0): If 1, then peptide sequence is randomized before starting PepBD design. If 0, then the sequence is not randomized
- SA_FLAG (default=0): If 1, then use simulated annealing to perform PepBD design. If 0, then use constant temperature.
- SA_INTERVAL (default=100): Number of PepBD steps to perform before reducing reference temperature for sequence and backbone changes
- SA_FRACTION (default=0.9): Factor by which reference temperature is multipled after each SA_INTERVAL
- SA_T_START (default=10): Starting reference system temperature for both sequence and backbone changes. Units of kcal/mol
- SA_T_END (default=0.5): Final reference system temperature for both sequence and backbone changes. Units of kcal/mol

## sys.pdb
Contains the peptide-receptor complex. Since PepBD_Plastics takes longer as the system size increases, typically receptor atoms that are 10 Angstroms from the peptide are removed, e.g. using VMD.
`sys.pdb` should be formatted such that
 - Only lines for atomic coordinates remain
 - The peptide appears first, followed by the receptor 
  
