# PepBD_Plastics
PepBD code used for designing plastic-binding peptides

### REQUIRED FILES FOR RUNNING PEPBD ###
The following files are required to run PepBD. Each is described in more detail below. All files should be placed in the directory for running PepBD, unless otherwise noted.
- PepBD_V.1.9.f90
- input.txt (see details below)
- a pdb file of the peptide-receptor (see details below)
- the lib folder, placed at the location ~/PepBD/
- An empty directory named "pdbfiles" in the same directory the program is being run

### REQUIRED PEPBD FILE - INPUT FILE ###
The name of this file must be "input.txt". An example input file is provided in the repository. The variables that can be specified in the input file are the following:
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

### REQUIRED PEPBD FILE - SYSTEM PDB FILE ###
Contains the peptide-receptor complex. Since PepBD runs slower as the system size increases, receptor atoms that are far from the peptide are removed from the pdb file. Common cutoff distances used are 8-10 Angstroms. 

After reducing the system size, prepare the pdb file so it can be understood by PepBD by doing the following
 - Remove any lines that do not have atomic coordinates, including TER and END
 - If needed, reorder the atoms so the peptide appears first in the file, then is followed by the receptor
    - If this is done, then fix the residue and atom numbering to match the new order. The tleap Amber of module is useful for this.

### REQUIRED PEPBD FILE - PepBD_V.1.7.f90 ###
Count the total number of residues in the final pdb file, then replace "gnum=0" in with "gnum={NUM_RESIDUES} in the PepBD_1.7.f90 file.

### COMPILING ###
To compile the program, the intel compiler needs to be used (the program does not run properly when compiled with gfortran!). Compiling appears to work (at least) with compiler versions 2017 to present. PepBD can be compiled via the command 

ifort -O2 -o main PepBD_V1.7.f90

Note that the program needs to be recompiled each time the system size changes, due to the "gnum" parameter changing. This has been fixed in more recent versions, but is retained here since this version of the program was used in the paper.

### RUNNING THE PROGRAM ###

Run the executable as normal. For example, if the executable is named "main", then it can be run from the terminal using
./main

The duration depends on the system size and the number of steps to be performed, but normally takes a day to a week to complete. Note that PepBD_V1.7.f90 is not parallelized and does not support GPU.

#### GENERATED DATA FILES ####
The outputs from PepBD are the following:
	-energyprofile.txt: score and energy components at each PepBD step
	-energydetails.txt: sequence, energy, and move type at each PepBD step
	-output.txt: Details of the PepBD run
	-"pdbfiles" directory, which contains pdbs periodically generated during PepBD design.
	-rmsd.txt: the rmsd relative to the starting state at each step of design. This is disabled in PepBD_V1.7.f90, so it outputs 0 at every step.
 
