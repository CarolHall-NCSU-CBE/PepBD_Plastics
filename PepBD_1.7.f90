module constant

    implicit none
! The Data type "groupdetails" is used to define an array to store the information on the residues and atoms.
	type groupdetails
		integer*8			:: cnum1, cnum2, cnum3
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*5			:: gtype
		double precision				:: coo1(20,3), coo2(60,3), coo3(20,3)
	end type

	type energyparameters
		integer*8			:: atomid1(20), atomid2(60), atomid3(20)
		double precision	:: charge1(20), epsion1(20), r1(20), rborn1(20), fs1(20), dielecons1(20)
		double precision				:: charge2(60), epsion2(60), r2(60), rborn2(60), fs2(60), dielecons2(60)
		double precision				:: charge3(20), epsion3(20), r3(20), rborn3(20), fs3(20), dielecons3(20)
	end type

	type databackup
		double precision				:: coo(3)
	end type

	type lib4aminoacid
		integer*8			:: cnum1, cnum2, cnum3
		integer*8			:: dihedralangle(35,6)
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*4			:: gtype
		double precision	:: coo1(20,3), coo2(60,3), coo3(20,3)
		integer*8			:: grade, rotanum
	end type

	type index4sidechain
		integer*8			:: class_No, member_No
	end type

	type conformer4sidechain
		double precision				:: member(15,3)
	end type

	type dihedralparameters
		integer*8 			:: iph(36), jph(36), kph(36), lph(36), multiply(36)
		double precision    			:: pk(36,4), pn(36,4), phase(36,4)
	end type

	type atomlink      ! The Data type "atomlink" is used to store the neighboring atoms for each atom.
		integer*8			:: linknum
		integer*8			:: linkindex(4)
	end type

	!integer*8					:: gnum    ! Total number of monomers in the system (peptide + receptor)
	!integer*8					:: atom_num
	!integer*8, dimension(:), allocatable	:: residuesite

	integer*8, parameter	:: gnum=0  ! Total number of monomers in the system (peptide + receptor)
	integer*8, parameter	:: atom_num=gnum*60 !stores info on all atoms in the system (factor of 60 arbitrarily chosen to be "big enough" to store everything, may want to fix this eventually)
	integer*8				:: residuesite(gnum) !stores atom indices for peptide residues (no receptor info stored)
	integer*8				:: nstep_start, nstep_terminal, idum, sitenum, recalcu_switch, rec_calculation_flag, flag4conformer
	integer*8				:: flavoredregion_number, rama_map(-179:180, -179:180) ! The defined varivables are for the ramachandran map.
	integer*8				:: fragmentnum !number of freely rotatable residues in the peptide
	integer*8				:: helix_num !number of helices in the peptdie
	integer*8, parameter	:: maximum_helix_num=5 !max number of helices in the peptide
	integer*8, parameter	:: number4peptide_candidates=6 !max number of conformation changes found by CONROT before it evaluats which of the 6 gives the lowest energy
	integer*8				:: helix_start(maximum_helix_num), helix_end(maximum_helix_num) !residues numbers for start/end of all helices
	integer*8				:: steps_before_conf_change = 50 !number of steps taken in each cycle before conformation changes are allowed
	double precision					:: ekt_seq = 1, ekt_backbone_ratio = 1,  ekt_backbone,  ekt_rigid_ratio =1,  ekt_rigid, ekt_scmf=0.6 !reference temperatures for sequence and conformation changes, and temperature for SCMF method. These can be changed in input file
	double precision					:: ph_value        ! The setting of PH
	double precision, parameter		:: lampda=0.5  ! The convergence coefficient
	double precision, parameter		:: surftens=0.0072  ! The surface tension
	double precision, parameter		:: offsetvalue=0.0
	double precision					:: weighting_factor=0.010  ! weighting factor is used to adjust the importance to the ligand within the score function
	double precision, parameter		:: dihedral_weighting_factor=0.50 !weights dihedral energy in score function
	double precision, parameter		:: vdw14_coeff=2.0  !1-4 VDW scale factor
	double precision, parameter		:: ele14_coeff=1.2  !1-4 ELE scale factor
	double precision					:: energy_min=10000 !stores minimum energy of peptide found during design (set to an arbitrary large positive number initially)
	double precision					:: rmsd_max    !The maximum restraints of rmsd
	double precision					:: rmsd_min    !The minimum restraints of rmsd
	double precision, parameter		:: pi=3.14159265358979323846264 !constant pi
	integer*8				:: PDB_gen_freq = 50 !PDB files generated after this many steps during design 
	integer*8				:: one=1, zero=0 !used in a couple locations to avoid integer size errors
	integer*8				:: rmsd_flag = 0 !tells program when rmsd should be calculated. 1 = calculate, 0 = don't calculate

    !Parameters for CONROT
    real                  :: delta_phi_min = -6.0, delta_phi_max = 6.0, delta_phi_step = 1.0

    !parameters for random translation/rotation
    real                  :: translate_max = 2.0, translate_min = 0.5 !min/max translation amount in single step
	real					:: rotate_max = pi !max rotation amount
	double precision					:: max_dev_from_center = 5 !number of Angstroms peptide's center of mass can deviate from original peptide center of mass during design
	real					:: start_pep_com(3) !stores starting center of mass coordinates of the peptide
	integer*8, parameter	:: num_top_layer = 100 !this many atoms are used to calculate the top of the polymer surface 
	integer*4				:: translate_z_flag = 0 !if 1, then uses receptor surface coordinates and max_dev_from_surface to determine peptide translations in z-direction. if 0, then use max_dev_from_center
	double precision					:: min_zdist_from_surface = 3, max_zdist_from_surface = 1000 !minimum and maximum allowed distance between peptide center of mass and surface (reqires surface normal to be in positive z-direction)	

	!probabilities governing peptide changes
	real					:: backbone_prob=0.5, backbone_switch  !Probability of doing a backbone change, and corresponding number for 
	real					:: seq_prob = 0.5, seq_switch !Probability of doing a sequence change.
	real					:: scmfswitch=0.8  ! probabililty of doing a residue mutation. (1-scmfswitch) is probability of exchanging two residues
	double precision					:: CONROT_two_term = 0.1, CONROT_one_term = 0.4
	real					:: rotation_switch = 0.5 !probability of doing rigid body rotation, vs. a rigid body translation

    !Parameters for fixing the number of tryptophans
    integer*8               :: trp_limit = 3, trp_count !default to only allowing 3 tryptophans in a peptide

    !Parameters for variable hydration
	integer*8 				:: min_hyd(6), max_hyd(6), cur_hyd(6), above_hyd(6), below_hyd(6), temp_hyd(6) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative

	!Variables for tracking mutation probability success
	integer*8				:: SEQ_SCMF_attempts=0, SEQ_POINT_attempts=0, CONF_N2C_attempts=0, CONF_N_attempts=0, CONF_C_attempts=0
	integer*8				:: CONF_Internal_attempts=0, rotate_attempts=0, translate_attempts=0
	integer*8, allocatable  :: SEQ_PerRes_attempts(:), SEQ_PerRes_success(:)
	integer*8				:: SEQ_SCMF_success=0, SEQ_POINT_success=0, CONF_N2C_success=0, CONF_N_success=0, CONF_C_success=0
	integer*8				:: CONF_Internal_success=0, rotate_success=0, translate_success=0

	!Parameters for restart method
	integer*8				:: num_restarts, steps_per_cycle, randomize_length=10
	integer*8				:: initial_randomize = 0 !determines whether peptide sequence is randomized before starting design; 1 = randomize, 0 = don't randomize
	integer*8 				:: steps_before_randomize=50 !how many steps will be taken with initial sequence before randomizing peptide, if initial_randomize == 0
	integer*8				:: smart_restart_flag = 0 !controls whether program actively tracks if a new cycle should be started
	integer*8				:: smart_restart_len = 100 !number of steps to average the scores over
	double precision					:: smart_restart_offset=5.0 !max offset between current score and best score to keep searching
    integer*8               :: smart_restart_wait_len = 5 !wait this many steps after randomizing sequence to start accumulating score (scores at beginning will be very high)
	integer*8				:: smart_restart_score_method = 1 !method for determining when to restart. 1 = compare average scores of intervals, 2 = compare best scores of intervals, 3 = compare average score of current interval to best score

	!parameters for simulated annealing --> currently using geometric series for temeprature schedule (T_{t+1}  = T_{t} * r, 0 < r < 1)
	integer*8				:: SA_flag = 0 !default to no simulated annealing
	integer*8				:: SA_interval = 100 !after SA_interval steps, reduce the temperature
	double precision					:: SA_fraction = 0.9 !value of r
	double precision					:: SA_T_start = 10, SA_T_end = 0.5 !default starting and end temperature
	integer*8				:: SA_num_cycles !number of intervals needed to go from SA_T_start to SA_T_end (will overshoot if necessary)
	integer*8 				:: SA_best_hyd(6), SA_above_hyd(6), SA_below_hyd(6), SA_trp_count !stores hydration properties of best peptide in current cycle

	!variables for per-residue binding energy
	double precision, allocatable		:: energy_per_res(:, :) !stores binding energy per residue
	double precision, allocatable		:: res_select_p_vector(:) !stores probability of selecting a residue to be mutated
	integer*8				:: decomp_mutate_flag = 0 !if this flag is 0, then select residues to mutate with likelihood based on their contribution to binding energy
	character*4				:: decomp_score_type = "INV"			

	!variables for non-polar solvation energy
	integer*8					:: ESURF_flag = 0 !determines if non-polar solvation energy is calculated
	integer*8					:: ESURF_count = 50	!number of top peptides to calculate ESURF for, default to 50 but can specified by user
	double precision, allocatable			:: top_pep_scores(:) !stores ordered scores for top ESURF_count peptides (order is lowest score at index 1, highest score at index N)
	character*15, allocatable 	:: ESURF_pep_descriptors(:) !stores prefix of pdb file for top peptides
	integer*8					:: ESURF_worst_index !stores index to worst scoring peptide in current top scores list
	double precision						:: ESURF_worst_score !stores score of worst scoring peptide in top_pep_scores

	!variables for neighbor list (currently not used)
    integer*8               :: neighbor_flag = 0 !if 1, then use neighbor cutoff. If 0, then don't use neighbor cutoff
	double precision					:: neighbor_cutoff = 1000.0 !beyond this distance, two atoms not considered neighbors and will not calculate pair-wise energy. By default, this value is large so all atoms are neighbors
    double precision					:: neighbor_cutoff2 !stores neighbor_cutoff squared, which is used frequently in calculations
    integer*8               :: neighbor_length = 200 !should be more than enough to store all neighbors for each residue in the peptide
    integer*8, allocatable  :: neighbor_list(:,:)   !will store index to neighbor residues for each residue in the ligand (NOTE: doing this on a per-residue basis) 
    integer*8, allocatable  :: neighbor_list_lens(:) !will store how many neighbors each ligand residue has
    real *8, allocatable    :: neighbor_rec_res_box_edges(:,:) !will store edges of each receptor's box edges to determine distance from each ligand residue' box edges

	!Tabu search parameters
	integer					:: TABU_flag = 0 !if 1, then tabu search used for sequence changes
	integer, allocatable	:: TABU_allowed_residues(:) !list of amino acid positions that CAN be mutated
	integer					:: TABU_ban_length = 0 !number of mutations that mutated residue cannot be changed. NOT EQUAL TO STEPS!! There are n sequence mutations per step, where n is the number of amino acids in the peptide
	integer, allocatable	:: TABU_ban_lengths(:) !stores list of how long tabu residues have been on the tabu list. list decremented, so residue can be mutated again once the length = 0
	integer					:: TABU_allowed_length !stores number of residues currently allowed by TABU search

	!other parameters
	double precision					:: rec_sgb=0.0, rec_snp=0.0 !stores energy info on the receptor 
	double precision					:: flavored_region(60000,2) ! It is for the ramachandran map.
	double precision                  :: zero_real = 0.0, one_real = 1.0	!zero and one required by some calls to subroutines
	integer					:: move_pep_over_surface_flag = 0 !if 1, then translate peptide in z-direction over surface before starting design. if 0, then do nothing
	double precision					:: target_pep_distance = 4.5 !if translating peptide over surface, this value is target distance between peptide center of mass and top of surface 
	double precision					:: pep_com(3), pol_com(3)
	real					:: surface_z, pep_min_z=1000000 !stores z coordinate of surface and minimum z coordinate of peptide

	character*30			:: filename1, filename2, receptor_name

type(groupdetails)      :: original_group(gnum)
type(lib4aminoacid)     :: aa_lib(20)

end module constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module utilities
	use constant
	contains
	
    recursive subroutine Find_Worst_Index(start, end, worst_score, worst_index,factor)
	implicit none
	double precision			:: worst_score, score1, score2
	integer*8		:: start, end, split, worst_index, index1, index2, factor
	
	if (start.ge.end) then
		worst_index = start
		worst_score = top_pep_scores(start)
	else
		split = int(start+(end-start)/2)
		call Find_Worst_Index(start, split, score1, index1, factor)
		call Find_Worst_Index(split+1, end, score2, index2, factor)
		if ((score1*factor).ge.(score2*factor)) then
			worst_index = index1
			worst_score = score1
		else
			worst_index = index2
			worst_score = score2
		endif
	endif
	return

	end subroutine Find_Worst_Index

	subroutine update_backup_coords(group, groupdata_backup)
	
		implicit none	
	type(groupdetails)			:: group(gnum)
	type(databackup)			:: groupdata_backup(gnum)
	integer*4					:: i,j 

	do i=1,sitenum
		do j=1,group(i)%cnum1
			if(group(i)%atype1(j)=="H".or.group(i)%atype1(j)=="H1") then
				groupdata_backup(i)%coo(:)=group(i)%coo1(j,:)
			endif
		enddo
	enddo

end subroutine update_backup_coords


end module utilities

module randomgenerator

	use constant
	use utilities

	contains
! Generate the random number
	subroutine ran_gen(ran2,flag_ran)
	implicit none
    integer*8			:: im1,im2,imm2,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv, start_flag
	integer*8				:: flag_ran
    double precision			:: ran2,am,eps,rnmx
    parameter       (im1=2147483563,im2=2147483399,am=1./im1,imm2=im1-1,ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211, &
					ir2=3791,ntab=32,ndiv=1+imm2/ntab,eps=1.2e-7,rnmx=1.-eps)
	integer*8			:: idum2,j,k,iv(ntab),iy
    save			iv,iy,idum2
	save			start_flag
    data            idum2/123456789/, iv/ntab*0/, iy/0/, start_flag/0/

	if(recalcu_switch.ne.0.and.start_flag.eq.0) then
		start_flag=1
		open(20, file="randomnumber.txt", status="old")
			read(20,*)
			read(20,*) idum
			read(20,*)
			read(20,*) idum2
			read(20,*)
			do j=1, ntab
				read(20,*) iv(j)
			enddo
			read(20,*)
			read(20,*) iy
		close(20)
	endif

	if(flag_ran==0) then
		if (idum.le.0) then
			idum=max(-idum,1)
			idum2=idum
			do 100 j=ntab+8,1,-1
				k=idum/iq1
				idum=ia1*(idum-k*iq1)-k*ir1
				if (idum.lt.0) idum=idum+im1
				if (j.le.ntab) iv(j)=idum
100			continue
			iy=iv(1)
		endif
		k=idum/iq1
		idum=ia1*(idum-k*iq1)-k*ir1
		if (idum.lt.0) idum=idum+im1
		k=idum2/iq2
		idum2=ia2*(idum2-k*iq2)-k*ir2
		if (idum2.lt.0) idum2=idum2+im2
		j=1+iy/ndiv
		iy=iv(j)-idum2
		iv(j)=idum
		if(iy.lt.1) iy=iy+imm2
		ran2=min(am*iy,rnmx)
	else
		open(20, file="randomnumber.txt", status="replace")
			write(20,*) "idum="
			write(20,*) idum
			write(20,*) "idum2="
			write(20,*) idum2
			write(20,*) "iv(ntab=32)="
			do j=1, ntab
				write(20,*) iv(j)
			enddo
			write(20,*) "iy="
			write(20,*) iy
		close(20)
	endif

    return
    end subroutine ran_gen

end module randomgenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input

    use constant

    contains
    subroutine inputfile
		implicit none
		integer*8             :: i=0, status, str_status, check_hyd
		character*40        :: name, val
	
		filename1 = 'temp'
		filename2 = 'temp'
		recalcu_switch = -1
		nstep_start = -1
		nstep_terminal = -1
		sitenum =-1
		receptor_name='temp'
		ph_value = -1
		rmsd_max = -1
		rmsd_min = -1
		min_hyd(1:6) = (/-1, -1, -1, -1, -1, -1/) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative
		max_hyd(1:6) = (/-1, -1, -1, -1, -1, -1/) !order is glycine, hydrophobic, positive, hydrophilic, other, and negative
	
		!read input file
		open(10, file="input.txt", status="old")
			do while(.true.)
				read(10, *, iostat=status) name, val
				if(status.ne.0) exit !end of file has been reached
	  
				if(name=="PDBFILE") then
					filename1 = val
				elseif(name=="RESTARTPDB") then
					filename2 = val
				elseif(name=="RECALCU_SWITCH") then
					read(val,*,iostat=str_status) recalcu_switch
				elseif(name=="STEP_START") then
					read(val,*,iostat=str_status) nstep_start
				elseif(name=="STEP_END") then
					read(val,*,iostat=str_status)  nstep_terminal
				elseif(name=="SITENUM") then
					read(val,*,iostat=str_status)  sitenum
				elseif(name=="RECEPTOR_NAME") then
					read(val,*,iostat=str_status) receptor_name
				elseif(name=="BACKBONE_CHANGE_PROB") then
					read(val,*,iostat=str_status) backbone_prob
				elseif(name=="SEQ_CHANGE_PROB") then
					read(val,*,iostat=str_status) seq_prob
				elseif(name=="SCMF_SWITCH") then
					read(val,*,iostat=str_status) scmfswitch
				elseif(name=="KT_SEQUENCE") then
					read(val,*,iostat=str_status) ekt_seq
				elseif(name=="KT_BACKBONE_RATIO") then
					read(val,*,iostat=str_status) ekt_backbone_ratio
				elseif(name=="KT_RIGID_RATIO") then
					read(val,*,iostat=str_status) ekt_rigid_ratio
				elseif(name=="PH") then
					read(val,*,iostat=str_status) ph_value
				elseif(name=="RANSEED") then
					read(val,*,iostat=str_status) idum
				elseif(name=="MAX_RMSD") then
					read(val,*,iostat=str_status) rmsd_max
				elseif(name=="MIN_RMSD") then
					read(val,*,iostat=str_status) rmsd_min
				elseif(name=="NUM_GLY_MIN") then
					read(val,*,iostat=str_status) min_hyd(1)
				elseif(name=="NUM_PHOBIC_MIN") then
					read(val,*,iostat=str_status) min_hyd(2)
				elseif(name=="NUM_POS_MIN") then
					read(val,*,iostat=str_status) min_hyd(3)
				elseif(name=="NUM_POLAR_MIN") then
					read(val,*,iostat=str_status) min_hyd(4)
				elseif(name=="NUM_OTHER_MIN") then
					read(val,*,iostat=str_status) min_hyd(5)
				elseif(name=="NUM_NEG_MIN") then
					read(val,*,iostat=str_status) min_hyd(6)
				elseif(name=="NUM_GLY_MAX") then
					read(val,*,iostat=str_status) max_hyd(1)
				elseif(name=="NUM_PHOBIC_MAX") then
					read(val,*,iostat=str_status) max_hyd(2)
				elseif(name=="NUM_POS_MAX") then
					read(val,*,iostat=str_status) max_hyd(3)
				elseif(name=="NUM_POLAR_MAX") then
					read(val,*,iostat=str_status) max_hyd(4)
				elseif(name=="NUM_OTHER_MAX") then
					read(val,*,iostat=str_status) max_hyd(5)
				elseif(name=="NUM_NEG_MAX") then
					read(val,*,iostat=str_status) max_hyd(6)
				elseif(name=="NUM_HELIX") then
					read(val,*,iostat=str_status) helix_num
					if(helix_num.gt.maximum_helix_num) then
						open(20, file="error.txt", access="append")
							write(20,*) "Error in input file!!"
							write(20,*) "helix_num=", helix_num, "is more than the maximum helix number", maximum_helix_num
							write(20,*) "Please adjust the maximum_helix_num in the code!"
						close(20)
						stop
					endif
				elseif(name=="HELIX_START") then
					read(val,*,iostat=str_status) helix_start(i)
				elseif(name=="HELIX_END") then
					read(val,*,iostat=str_status) helix_end(i)
					i = i + 1
				elseif(name=="TRP_LIMIT") then
					read(val,*,iostat=str_status) trp_limit
				elseif(name=="STEP_BEFORE_CONF_CHANGE") then
					read(val,*,iostat=str_status) steps_before_conf_change	
				elseif(name=="NUM_RESTARTS") then
					read(val,*,iostat=str_status) num_restarts
				elseif(name=="STEPS_PER_CYCLE") then
					read(val,*,iostat=str_status) steps_per_cycle
				elseif(name=="RANDOMIZE_LENGTH") then
					read(val,*,iostat=str_status) randomize_length
				elseif(name=="PDB_GEN_FREQ") then
					read(val,*,iostat=str_status) PDB_gen_freq
				elseif(name=="PEP_WEIGHT_FACTOR") then
					read(val,*,iostat=str_status) weighting_factor
				elseif(name=="CONROT_TWO_TERM") then
					read(val,*,iostat=str_status) CONROT_two_term
				elseif(name=="CONROT_ONE_TERM") then
					read(val,*,iostat=str_status) CONROT_one_term
				elseif(name=="INITIAL_RANDOMIZE") then
					read(val,*,iostat=str_status) initial_randomize	
				elseif(name=="STEPS_BEFORE_RANDOMIZE") then
					read(val,*,iostat=str_status) steps_before_randomize	
				elseif(name=="SMART_RESTART_FLAG") then
					read(val,*,iostat=str_status) smart_restart_flag
				elseif(name=="SMART_RESTART_LEN") then
					read(val,*,iostat=str_status) smart_restart_len	
				elseif(name=="SMART_RESTART_OFFSET") then
					read(val,*,iostat=str_status) smart_restart_offset
				elseif(name=="SMART_RESTART_WAIT_LEN") then
					read(val,*,iostat=str_status) smart_restart_wait_len
				elseif(name=="SMART_RESTART_SCORE_METHOD") then
					read(val,*,iostat=str_status) smart_restart_score_method
				elseif(name=="SA_FLAG") then
					read(val,*,iostat=str_status) SA_flag
				elseif(name=="SA_INTERVAL") then
					read(val,*,iostat=str_status) SA_interval	
				elseif(name=="SA_FRACTION") then
					read(val,*,iostat=str_status) SA_fraction
				elseif(name=="SA_T_START") then
					read(val,*,iostat=str_status) SA_T_start
				elseif(name=="SA_T_END") then
					read(val,*,iostat=str_status) SA_T_end
				elseif(name=="DECOMP_MUTATE_FLAG") then
					read(val,*,iostat=str_status) decomp_mutate_flag
				elseif(name=="DECOMP_SCORE_TYPE") then
					read(val,*,iostat=str_status) decomp_score_type
				elseif(name=="NEIGHBOR_FLAG") then
					read(val,*,iostat=str_status) neighbor_flag
				elseif(name=="NEIGHBOR_CUTOFF") then
					read(val,*,iostat=str_status) neighbor_cutoff
				elseif(name=="ROTATION_SWITCH") then
					read(val,*,iostat=str_status) rotation_switch
				elseif(name=="TRANSLATE_Z_FLAG") then
					read(val,*,iostat=str_status) translate_z_flag
				elseif(name=="MIN_ZDIST_FROM_SURFACE") then
					read(val,*,iostat=str_status) min_zdist_from_surface
				elseif(name=="MAX_ZDIST_FROM_SURFACE") then
					read(val,*,iostat=str_status) max_zdist_from_surface
				elseif(name=="TABU_FLAG") then
					read(val,*,iostat=str_status) TABU_flag
				elseif(name=="TABU_BAN_LENGTH") then
					read(val,*,iostat=str_status) TABU_ban_length
				elseif(name=="MOVE_PEP_OVER_SURFACE") then
					read(val,*,iostat=str_status) move_pep_over_surface_flag
				elseif(name=="PEP_OVER_SURFACE_Z") then
					read(val,*,iostat=str_status) target_pep_distance
				else
					open(20, file="error.txt", access="append")
					write(20,*) "Error in input file!!"
					write(20,*) "Unrecognized parameter ", name, " in input file, terminating program."
					close(20)
					stop
				endif
		end do
	
		!check the input file is okay 
		open(20, file="error.txt", access="append")
	
		!check enough helix start/stop positions were given for each helix
		if (i.ne.helix_num) then
			write(20,*) "Error in input file!!"
			write(20,*) "Incorrect number of helix start/stop positions given, terminating program."
			write(20,*) "i is ", i, " while there are ", helix_num, "helices"
			close(20)
			stop
		endif
	  
		!check necessary inputs were provided
		if (filename1 == 'temp'.or.filename2 == 'temp'.or.recalcu_switch == -1.or.nstep_start == -1.or. &
			sitenum == -1.or.receptor_name == 'temp'.or. &
			ph_value == -1.or.rmsd_max == -1.or.rmsd_min == -1.or.min_hyd(1)== -1.or.&
			min_hyd(2) == -1.or.min_hyd(3) == -1.or.min_hyd(4) == -1.or.min_hyd(5) == -1.or.min_hyd(6) == -1.or. &
			max_hyd(1) == -1.or.max_hyd(2)== -1.or.max_hyd(3) == -1.or.max_hyd(4) == -1.or.max_hyd(5)== -1.or.max_hyd(6)==-1) then
			write(20,*) "Error in input file!!"
			write(20,*) "Not all required inputs were provided, please check input file. Terminating program."
			close(20)
			stop
		end if
	
		!check the restart parameters make sense
		if (num_restarts.eq.-1) then 
			if (nstep_terminal.eq.-1) then
				write(20,*) "Error in input file!!"
				write(20,*) "Please either provide restart parameters (NUM_RESTARTS and STEPS_PER_CYCLE) &
							or length of single run (STEP_END). Ending program"
				close(20)
				stop
			else
				num_restarts=1 
				steps_per_cycle = nstep_terminal - nstep_start
			endif
		endif
	
		if (num_restarts.gt.1.and.steps_per_cycle.eq.-1) then
			write(20,*) "Error in input file!!"
			write(20,*) "Trying to use restart method, but number of steps before restarting not specified!"
			write(20,*) "Add STEPS_PER_CYCLE parameter to input.txt file"
			close(20)
			stop
		endif
	
		!check the hydration properties can be used
		check_hyd = max_hyd(1) + max_hyd(2) + max_hyd(3) + max_hyd(4) + max_hyd(5) + max_hyd(6) 
		if (check_hyd < sitenum) then
			write(20,*) "Error in input file!!"
			write(20,*) "The sum of the maximum hydration properties, ", check_hyd, "is less than the &
			 total number of peptide amino acids, ", sitenum, " so design is not possible!"
			write(20,*) "Fix the maximum hydration values so it at least equals the number of amino acids in the peptide"
			close(20)
			stop
		endif
	
		!check the move probabilities are usable
		if (CONROT_two_term.ge.CONROT_one_term) then
			write(20,*) "Error in input file!!"
			write(20,*) "CONROT_two_term probability is larger than CONROT_one_term! Please fix input file so CONROT_one_term is larger"
			close(20)
			stop
		endif
	
		if (backbone_prob + seq_prob > 1.00001) then
			write(20,*) "Error in input file!!"
			write(20,*) "Sum of BACKBONE_PROB and SEQ_PROB cannot be greater than 1! Please fix the input file"
			close(20)
			stop
		elseif(backbone_prob.lt.0.0.or.seq_prob.lt.0.0) then
			write(20,*) "Error in input file!!"
			write(20,*) "Neither BACKBONE_PROB nor SEQ_PROB can be less than 0! Please fix the input file"
			close(20)
			stop
		else
			seq_switch = 1.0-seq_prob
			backbone_switch = seq_switch - backbone_prob
		endif
	
		!check the smart restart inputs
		if (smart_restart_flag.eq.1.and.smart_restart_len.gt.steps_per_cycle) then
			write(20,*) "Error in input file!!"
			write(20,*) "value of SMART_RESTART_LEN cannot be greater than STEPS_PER_CYCLE! Please update input file and rerun"
			close(20)
			stop
		endif
	
		if (smart_restart_score_method.lt.1.or.smart_restart_score_method.gt.3) then
			write(20,*) "Error in input file!!"
			write(20,*) "value of SMART_RESTART_SCORE_METHOD not usable! This value can only be 1, 2, or 3. Please update input file and rerun"
			close(20)
			stop
		endif
	
		!check the simulated annealing inputs
		if (SA_flag.eq.1) then
			if (SA_T_start.le.SA_T_end) then  
				write(20,*) "Error in input file!!"      
				write(20,*) "Starting temperature for simulated annealing is smaller than end temperature! Please update input file and rerun"
				close(20)
				stop
			endif
			if (smart_restart_flag.eq.1) then
				write(20,*) "Error in input file!!"
				write(20,*) "Currently do not allow both simulated annealing and random restart. Please update input file to remove one of these options and rerun"
				close(20)
				stop
			end if
	
			!calculate number of intervals needed for annealing
			SA_num_cycles = int((log(SA_T_end) - log(SA_T_start))/log(SA_fraction)) + 1
			steps_per_cycle = SA_interval
		endif
	
		!other miscellaneous checks
		if (translate_z_flag.ne.0.and.translate_z_flag.ne.1) then
			write(20,*) "Error in input file!!"
			write(20,*) "TRANSLATE_Z_FLAG can only be 0 or 1. Please update input file and rerun."
			close(20)
			stop
		endif
	
		if (min_zdist_from_surface.lt.0) then
			write(20,*) "Error in input file!!"
			write(20,*) "MIN_ZDIST_FROM_SURFACE should not be less than 0, as this causes steric overlaps and results in many translations being rejected."
			write(20,*) "Please update input file to increase value above 0 and rerun."
			close(20)
			stop
		endif
	
		if (ph_value.lt.6.0.or.ph_value.gt.8.3) then
			write(20,*) "Error in input file!!"
			write(20,*) "Input pH is outside the range of 6.0 - 8.3, and the program may not work in this region."
			write(20,*) "This will be fixed in the future. Terminating program."
		endif
	
		if (decomp_mutate_flag.eq.1.and.(decomp_score_type.ne."REG".and.decomp_score_type.ne."INV")) then
			write(20,*) "Error in input file!!"
			write(20,*) "When using per-residue energy composition for mutations, only valid options for DECOMP_SCORE_TYPE are INV or REG."
			write(20,*) "Input value of ", decomp_score_type, " not allowed!"
			write(20,*) "Please update input file and rerun."
			close(20)
			stop
		endif
	
		if (TABU_flag.eq.1.and.decomp_mutate_flag.eq.1) then
			write(20,*) "Error in input file!!"
			write(20,*) "Currently cannot use both TABU search and select residues for mutation based on energy decomposition, as residues will not be selected properly"
			write(20,*) "Please update input file to use only one of these options and rerun."
			close(20)
			stop
		endif
	
		close(20)
	
		!neighbor cutoff terms
		neighbor_cutoff2 = neighbor_cutoff * neighbor_cutoff
	
			!setup reference temperatures
		ekt_backbone = ekt_seq * ekt_backbone_ratio
		ekt_rigid = ekt_seq * ekt_rigid_ratio
	
		!fix amino acid names based on pH of design
		!call Get_AA_List()
	
		!setup the TABU lists
		if (TABU_flag.eq.1) then
			call initialize_TABU(1)
		endif
	 
		return
	end subroutine inputfile

	subroutine initialize_TABU(flag)
		implicit none
		integer    :: i,flag
	
		if (flag.eq.1) then !only need to allocate memory the first time we initialize TABU list
			allocate(TABU_ban_lengths(sitenum))
			allocate(TABU_allowed_residues(sitenum))
		endif

		TABU_ban_lengths = 0
		do i=1,sitenum
			TABU_allowed_residues(i) = i
		enddo
	
		TABU_allowed_length = sitenum
	end subroutine initialize_TABU
	
	!subroutine allocate_memory(filename1, group, groupdata_backup, gnum, atom_num, W_numex, W_numex4, W_inb, W_inb4, group_para)
    !determine how many residues are in the system, then allocate memory
	!implicit none
	!integer*8,intent(OUT) 		:: gnum, atom_num
    !integer         :: status=0, anum
    !character*20	:: filename1
    !double precision		    :: x, y, z
	!character*4		:: char, atype, name
    !type(groupdetails), dimension(:), allocatable	:: group
	!type(databackup), dimension(:), allocatable		:: groupdata_backup
    !integer*8, dimension(:), allocatable :: W_numex, W_numex4
	!integer*8, dimension(:,:), allocatable :: W_inb, W_inb4
	!type(energyparameters), dimension(:), allocatable :: group_para

    !open(10, file=filename1)
    !do while(status.eq.0) 
    !    read(10, *, iostat=status) char, anum, atype, name, anum, x, y, z
    !enddo
    !close(10)

    !allocate(original_group(gnum))
    !allocate(group(gnum))
    !allocate(groupdata_backup(gnum))
    !allocate(W_numex(atom_num)); allocate(W_numex4(atom_num))
	!allocate( W_inb(atom_num,20)); allocate(W_inb4(atom_num, 60))
	!allocate(group_para(gnum))

	!end subroutine allocate_memory
  
end module input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdbfile

	use	constant
	use utilities

	contains

	subroutine readpdb_for_ESURF(group, gnum, fname)
	implicit none
	integer*8					:: anum, status=1, flag1, gnum
	double precision						:: x, y, z
	character*4					:: char, atype, name
	character*30				:: fname
	type(groupdetails)			:: group(gnum)

	group%cnum1=0
	group%cnum2=0
	group%cnum3=0
	flag1=0
	open(10, file=trim(fname))
	do while(.true.)
		read(10, *, iostat=status) char, anum, atype, name, anum, x, y, z
		if(status.lt.0) exit
		if (status.eq.0) then !when status > 0, then we are skipping over line that is the "TER" separating peptide from receptor
			if(anum.le.sitenum) then
				if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
				.or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
				.or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
				.or.name=="TYR".or.name=="VAL") then
					if (anum.eq.1) then
						name = "N" // name
					elseif (anum.eq.sitenum) then
						name = "C" // name
					endif
					group(anum)%gtype=name
					if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
						.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
						group(anum)%cnum1=group(anum)%cnum1+1
						group(anum)%atype1(group(anum)%cnum1)=atype
						group(anum)%coo1(group(anum)%cnum1,1)=x
						group(anum)%coo1(group(anum)%cnum1,2)=y
						group(anum)%coo1(group(anum)%cnum1,3)=z
					elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
						group(anum)%cnum3=group(anum)%cnum3+1
						group(anum)%atype3(group(anum)%cnum3)=atype
						group(anum)%coo3(group(anum)%cnum3,1)=x
						group(anum)%coo3(group(anum)%cnum3,2)=y
						group(anum)%coo3(group(anum)%cnum3,3)=z
					else
						group(anum)%cnum2=group(anum)%cnum2+1
						group(anum)%atype2(group(anum)%cnum2)=atype
						group(anum)%coo2(group(anum)%cnum2,1)=x
						group(anum)%coo2(group(anum)%cnum2,2)=y
						group(anum)%coo2(group(anum)%cnum2,3)=z
					endif
				endif
			else
				if(receptor_name=="rna") then
					if(name=="RC".or.name=="RA".or.name=="RU".or.name=="RG".or.name=="STA".or.name=="SMU".or.name=="1MA".or.name=="PSU"  &
						.or.name=="5MU".or.name=="RC5".or.name=="RA5".or.name=="RU5".or.name=="RG5".or.name=="RC3".or.name=="RA3" &
						.or.name=="RU3".or.name=="RG3".or.name=="2SU".or.name=="6TA") then
						group(anum)%gtype=name
						if(atype=="HO5'".or.atype=="H5T".or.atype=="P".or.atype=="O1P".or.atype=="OP1".or.atype=="O2P".or.atype=="OP2" &
						.or.atype=="O5'".or.atype=="C5'".or.atype=="1H5'".or.atype=="H5'".or.atype=="2H5'".or.atype=="H5''".or. &
						atype=="C4'".or.atype=="H4'".or.atype=="O4'".or.atype=="C1'".or.atype=="H1'".or.atype=="H5'1".or.atype=="H5'2") then
							if(atype=="H5'1") then
								group(anum)%cnum1=group(anum)%cnum1+1
								group(anum)%atype1(group(anum)%cnum1)="1H5'"
								group(anum)%coo1(group(anum)%cnum1,1)=x
								group(anum)%coo1(group(anum)%cnum1,2)=y
								group(anum)%coo1(group(anum)%cnum1,3)=z
							elseif(atype=="H5'2") then
								group(anum)%cnum1=group(anum)%cnum1+1
								group(anum)%atype1(group(anum)%cnum1)="2H5'"
								group(anum)%coo1(group(anum)%cnum1,1)=x
								group(anum)%coo1(group(anum)%cnum1,2)=y
								group(anum)%coo1(group(anum)%cnum1,3)=z
							else
								group(anum)%cnum1=group(anum)%cnum1+1
								group(anum)%atype1(group(anum)%cnum1)=atype
								group(anum)%coo1(group(anum)%cnum1,1)=x
								group(anum)%coo1(group(anum)%cnum1,2)=y
								group(anum)%coo1(group(anum)%cnum1,3)=z
							endif
						elseif(atype=="C3'".or.atype=="H3'".or.atype=="C2'".or.atype=="1H2'".or.atype=="H2'".or.atype=="O2'" &
							.or.atype=="2HO'".or.atype=="HO2'".or.atype=="H3T".or.atype=="HO3'".or.atype=="O3'".or.atype=="H2'1".or.atype=="HO'2") then
							if(atype=="H2'1") then
								group(anum)%cnum3=group(anum)%cnum3+1
								group(anum)%atype3(group(anum)%cnum3)="1H2'"
								group(anum)%coo3(group(anum)%cnum3,1)=x
								group(anum)%coo3(group(anum)%cnum3,2)=y
								group(anum)%coo3(group(anum)%cnum3,3)=z
							elseif(atype=="HO'2") then
								group(anum)%cnum3=group(anum)%cnum3+1
								group(anum)%atype3(group(anum)%cnum3)="2HO'"
								group(anum)%coo3(group(anum)%cnum3,1)=x
								group(anum)%coo3(group(anum)%cnum3,2)=y
								group(anum)%coo3(group(anum)%cnum3,3)=z
							else
								group(anum)%cnum3=group(anum)%cnum3+1
								group(anum)%atype3(group(anum)%cnum3)=atype
								group(anum)%coo3(group(anum)%cnum3,1)=x
								group(anum)%coo3(group(anum)%cnum3,2)=y
								group(anum)%coo3(group(anum)%cnum3,3)=z
							endif
						else
							group(anum)%cnum2=group(anum)%cnum2+1
							group(anum)%atype2(group(anum)%cnum2)=atype
							group(anum)%coo2(group(anum)%cnum2,1)=x
							group(anum)%coo2(group(anum)%cnum2,2)=y
							group(anum)%coo2(group(anum)%cnum2,3)=z
						endif
					else
						open(20, file="error.txt", access="append")
							write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
							write(20,*) "Please check whether the group type is right or not!"
						close(20)
						stop
					endif

				elseif(receptor_name=="peptide") then
					if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
						.or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
						.or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
						.or.name=="TYR".or.name=="VAL".or.name=="NALA".or.name=="NARG".or.name=="NASN".or.name=="NASP" &
						.or.name=="NCYS".or.name=="NGLN".or.name=="NGLU".or.name=="NGLY".or.name=="NHIE".or.name=="NILE" &
						.or.name=="NLEU".or.name=="NLYS".or.name=="NMET".or.name=="NPHE".or.name=="NPRO".or.name=="NSER" &
						.or.name=="NTHR".or.name=="NTRP".or.name=="NTYR".or.name=="NVAL".or.name=="CALA".or.name=="CARG" &
						.or.name=="CASN".or.name=="CASP".or.name=="CCYS".or.name=="CGLN".or.name=="CGLU".or.name=="CGLY" &
						.or.name=="CHIE".or.name=="CILE".or.name=="CLEU".or.name=="CLYS".or.name=="CMET".or.name=="CPHE" &
						.or.name=="CPRO".or.name=="CSER".or.name=="CTHR".or.name=="CTRP".or.name=="CTYR".or.name=="CVAL" &
						.or.name=="TYX".or.name=="ARN".or.name=="HIP".or.name=="CYT".or.name=="LYN".or.name=="GLH".or.name=="ASH" &
						.or.name=="NTYX".or.name=="NARN".or.name=="NHIP".or.name=="NCYT".or.name=="NLYN".or.name=="NGLH".or.name=="NASH" &
						.or.name=="CTYX".or.name=="CARN".or.name=="CHIP".or.name=="CCYT".or.name=="CLYN".or.name=="CGLH".or.name=="CASH") then
						group(anum)%gtype=name
						if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
							.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
							group(anum)%cnum1=group(anum)%cnum1+1
							group(anum)%atype1(group(anum)%cnum1)=atype
							group(anum)%coo1(group(anum)%cnum1,1)=x
							group(anum)%coo1(group(anum)%cnum1,2)=y
							group(anum)%coo1(group(anum)%cnum1,3)=z
						elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
							group(anum)%cnum3=group(anum)%cnum3+1
							group(anum)%atype3(group(anum)%cnum3)=atype
							group(anum)%coo3(group(anum)%cnum3,1)=x
							group(anum)%coo3(group(anum)%cnum3,2)=y
							group(anum)%coo3(group(anum)%cnum3,3)=z
						else
							group(anum)%cnum2=group(anum)%cnum2+1
							group(anum)%atype2(group(anum)%cnum2)=atype
							group(anum)%coo2(group(anum)%cnum2,1)=x
							group(anum)%coo2(group(anum)%cnum2,2)=y
							group(anum)%coo2(group(anum)%cnum2,3)=z
						endif
					else
						open(20, file="error.txt", access="append")
							write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
							write(20,*) "Please check whether the group type is right or not!"
						close(20)
						stop
					endif

				elseif(receptor_name=="molecule") then
					group(anum)%gtype=name
					group(anum)%cnum2=group(anum)%cnum2+1
					group(anum)%atype2(group(anum)%cnum2)=atype
					group(anum)%coo2(group(anum)%cnum2,1)=x
					group(anum)%coo2(group(anum)%cnum2,2)=y
					group(anum)%coo2(group(anum)%cnum2,3)=z
				endif
			endif
		endif
	enddo

	endsubroutine readpdb_for_ESURF
! Read an initial file.
	subroutine readpdb(group, groupdata_backup, gnum)
	implicit none
	integer*8						:: anum, status=1, flag, flag1, i, gnum
	double precision						:: x, y, z
	double precision						:: pol_top_z(num_top_layer), min_top_z 
	integer*8					:: min_top_z_index
	character*4					:: char, atype, name
	type(groupdetails)			:: group(gnum)
	type(databackup)			:: groupdata_backup(gnum)
	
	group%cnum1=0
	group%cnum2=0
	group%cnum3=0
	fragmentnum=0
	flag=0
	flag1=0
	pol_top_z = -1000
	min_top_z = -1000
	min_top_z_index = 1

	!first determine how many residues are in the system
	!open(10, file=filename1)
	!do while(status.ne.0)
	!		read(10, *, iostat=status) char, anum, atype, name, anum, x, y, z
		
	!enddo
	!gnum = anum

	!call allocate(gnum)

	flag1=0
	open(10, file=trim(filename1))
	do while(.true.)
		read(10, *, iostat=status) char, anum, atype, name, anum, x, y, z
		if(status.ne.0) exit
		if (anum.gt.gnum) then
			open(20, file='error.txt', access='append')
			write(20,*) "The number of residues in the pdb file is larger than gnum! Please update gnum in the code and rerun"
			stop
			close(20)
		endif
		if(anum.le.sitenum) then
			if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
			   .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
			   .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
			   .or.name=="TYR".or.name=="VAL") then
				if (anum.eq.1) then
					name = "N" // name
				elseif (anum.eq.sitenum) then
					name = "C" // name
				endif
				group(anum)%gtype=name
				if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
					.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
					group(anum)%cnum1=group(anum)%cnum1+1
					group(anum)%atype1(group(anum)%cnum1)=atype
					group(anum)%coo1(group(anum)%cnum1,1)=x
					group(anum)%coo1(group(anum)%cnum1,2)=y
					group(anum)%coo1(group(anum)%cnum1,3)=z
					if(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
						if(recalcu_switch==0) flag1=1
					else
						if(atype=="H".or.atype=="H1") then
							groupdata_backup(anum)%coo(1)=x
							groupdata_backup(anum)%coo(2)=y
							groupdata_backup(anum)%coo(3)=z
						endif
					endif
				elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
					group(anum)%cnum3=group(anum)%cnum3+1
					group(anum)%atype3(group(anum)%cnum3)=atype
					group(anum)%coo3(group(anum)%cnum3,1)=x
					group(anum)%coo3(group(anum)%cnum3,2)=y
					group(anum)%coo3(group(anum)%cnum3,3)=z
				else
					group(anum)%cnum2=group(anum)%cnum2+1
					group(anum)%atype2(group(anum)%cnum2)=atype
					group(anum)%coo2(group(anum)%cnum2,1)=x
					group(anum)%coo2(group(anum)%cnum2,2)=y
					group(anum)%coo2(group(anum)%cnum2,3)=z
				endif

				if(anum.ne.flag) then
					flag=anum
					do i=1, helix_num
						if(anum.ge.helix_start(i).and.anum.le.helix_end(i)) goto 10
					enddo
					fragmentnum=fragmentnum+1
					residuesite(fragmentnum)=anum
10					continue
				endif
			else
				open(20, file="error.txt", access="append")
					write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
					write(20,*) "Please check whether the group type is right or not!"
				close(20)
				stop
			endif

		else
			if(receptor_name=="rna") then
				if(name=="RC".or.name=="RA".or.name=="RU".or.name=="RG".or.name=="STA".or.name=="SMU".or.name=="1MA".or.name=="PSU"  &
		       		.or.name=="5MU".or.name=="RC5".or.name=="RA5".or.name=="RU5".or.name=="RG5".or.name=="RC3".or.name=="RA3" &
					.or.name=="RU3".or.name=="RG3".or.name=="2SU".or.name=="6TA") then
					group(anum)%gtype=name
					if(atype=="HO5'".or.atype=="H5T".or.atype=="P".or.atype=="O1P".or.atype=="OP1".or.atype=="O2P".or.atype=="OP2" &
					.or.atype=="O5'".or.atype=="C5'".or.atype=="1H5'".or.atype=="H5'".or.atype=="2H5'".or.atype=="H5''".or. &
					atype=="C4'".or.atype=="H4'".or.atype=="O4'".or.atype=="C1'".or.atype=="H1'".or.atype=="H5'1".or.atype=="H5'2") then
						if(atype=="H5'1") then
							group(anum)%cnum1=group(anum)%cnum1+1
							group(anum)%atype1(group(anum)%cnum1)="1H5'"
							group(anum)%coo1(group(anum)%cnum1,1)=x
							group(anum)%coo1(group(anum)%cnum1,2)=y
							group(anum)%coo1(group(anum)%cnum1,3)=z
						elseif(atype=="H5'2") then
							group(anum)%cnum1=group(anum)%cnum1+1
							group(anum)%atype1(group(anum)%cnum1)="2H5'"
							group(anum)%coo1(group(anum)%cnum1,1)=x
							group(anum)%coo1(group(anum)%cnum1,2)=y
							group(anum)%coo1(group(anum)%cnum1,3)=z
						else
							group(anum)%cnum1=group(anum)%cnum1+1
							group(anum)%atype1(group(anum)%cnum1)=atype
							group(anum)%coo1(group(anum)%cnum1,1)=x
							group(anum)%coo1(group(anum)%cnum1,2)=y
							group(anum)%coo1(group(anum)%cnum1,3)=z
						endif
					elseif(atype=="C3'".or.atype=="H3'".or.atype=="C2'".or.atype=="1H2'".or.atype=="H2'".or.atype=="O2'" &
						.or.atype=="2HO'".or.atype=="HO2'".or.atype=="H3T".or.atype=="HO3'".or.atype=="O3'".or.atype=="H2'1".or.atype=="HO'2") then
						if(atype=="H2'1") then
							group(anum)%cnum3=group(anum)%cnum3+1
							group(anum)%atype3(group(anum)%cnum3)="1H2'"
							group(anum)%coo3(group(anum)%cnum3,1)=x
							group(anum)%coo3(group(anum)%cnum3,2)=y
							group(anum)%coo3(group(anum)%cnum3,3)=z
						elseif(atype=="HO'2") then
							group(anum)%cnum3=group(anum)%cnum3+1
							group(anum)%atype3(group(anum)%cnum3)="2HO'"
							group(anum)%coo3(group(anum)%cnum3,1)=x
							group(anum)%coo3(group(anum)%cnum3,2)=y
							group(anum)%coo3(group(anum)%cnum3,3)=z
						else
							group(anum)%cnum3=group(anum)%cnum3+1
							group(anum)%atype3(group(anum)%cnum3)=atype
							group(anum)%coo3(group(anum)%cnum3,1)=x
							group(anum)%coo3(group(anum)%cnum3,2)=y
							group(anum)%coo3(group(anum)%cnum3,3)=z
						endif
					else
						group(anum)%cnum2=group(anum)%cnum2+1
						group(anum)%atype2(group(anum)%cnum2)=atype
						group(anum)%coo2(group(anum)%cnum2,1)=x
						group(anum)%coo2(group(anum)%cnum2,2)=y
						group(anum)%coo2(group(anum)%cnum2,3)=z
					endif
				else
					open(20, file="error.txt", access="append")
						write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
						write(20,*) "Please check whether the group type is right or not!"
					close(20)
					stop
				endif

			elseif(receptor_name=="peptide") then
				if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
					.or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
					.or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
					.or.name=="TYR".or.name=="VAL".or.name=="NALA".or.name=="NARG".or.name=="NASN".or.name=="NASP" &
					.or.name=="NCYS".or.name=="NGLN".or.name=="NGLU".or.name=="NGLY".or.name=="NHIE".or.name=="NILE" &
					.or.name=="NLEU".or.name=="NLYS".or.name=="NMET".or.name=="NPHE".or.name=="NPRO".or.name=="NSER" &
					.or.name=="NTHR".or.name=="NTRP".or.name=="NTYR".or.name=="NVAL".or.name=="CALA".or.name=="CARG" &
					.or.name=="CASN".or.name=="CASP".or.name=="CCYS".or.name=="CGLN".or.name=="CGLU".or.name=="CGLY" &
					.or.name=="CHIE".or.name=="CILE".or.name=="CLEU".or.name=="CLYS".or.name=="CMET".or.name=="CPHE" &
					.or.name=="CPRO".or.name=="CSER".or.name=="CTHR".or.name=="CTRP".or.name=="CTYR".or.name=="CVAL" &
					.or.name=="TYX".or.name=="ARN".or.name=="HIP".or.name=="CYT".or.name=="LYN".or.name=="GLH".or.name=="ASH" &
					.or.name=="NTYX".or.name=="NARN".or.name=="NHIP".or.name=="NCYT".or.name=="NLYN".or.name=="NGLH".or.name=="NASH" &
					.or.name=="CTYX".or.name=="CARN".or.name=="CHIP".or.name=="CCYT".or.name=="CLYN".or.name=="CGLH".or.name=="CASH") then
					group(anum)%gtype=name
					if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
						.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
						group(anum)%cnum1=group(anum)%cnum1+1
						group(anum)%atype1(group(anum)%cnum1)=atype
						group(anum)%coo1(group(anum)%cnum1,1)=x
						group(anum)%coo1(group(anum)%cnum1,2)=y
						group(anum)%coo1(group(anum)%cnum1,3)=z
					elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
						group(anum)%cnum3=group(anum)%cnum3+1
						group(anum)%atype3(group(anum)%cnum3)=atype
						group(anum)%coo3(group(anum)%cnum3,1)=x
						group(anum)%coo3(group(anum)%cnum3,2)=y
						group(anum)%coo3(group(anum)%cnum3,3)=z
					else
						group(anum)%cnum2=group(anum)%cnum2+1
						group(anum)%atype2(group(anum)%cnum2)=atype
						group(anum)%coo2(group(anum)%cnum2,1)=x
						group(anum)%coo2(group(anum)%cnum2,2)=y
						group(anum)%coo2(group(anum)%cnum2,3)=z
					endif
				else
					open(20, file="error.txt", access="append")
						write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
						write(20,*) "Please check whether the group type is right or not!"
					close(20)
					stop
				endif

			elseif(receptor_name=="molecule") then
				group(anum)%gtype=name
				group(anum)%cnum2=group(anum)%cnum2+1
				group(anum)%atype2(group(anum)%cnum2)=atype
				group(anum)%coo2(group(anum)%cnum2,1)=x
				group(anum)%coo2(group(anum)%cnum2,2)=y
				group(anum)%coo2(group(anum)%cnum2,3)=z

				!check if z should be added to the list of polymer surface atoms
				if (z.gt.min_top_z) then
					top_pep_scores(min_top_z_index) = z
					call Find_Worst_Index(one, num_top_layer,min_top_z,min_top_z_index, -1*one)
				endif
			endif
		endif
	enddo
	close(10)

	!calculate z-coordinate of top of polymer surface
	surface_z = sum(top_pep_scores)/num_top_layer

	if(anum.ne.gnum) then
		open(20, file="error.txt", access="append")
			write(20,*) "gnum does not equal total number of residues in the system!"
			write(20,*) "Update gnum in the program, recompile, and rerun"
		close(20)
		stop
	endif
	if(fragmentnum.lt.3) then
		open(20, file="error.txt", access="append")
			write(20,*) "The code works only for the situation there are at leat 3 amino &
			acids as the pivots to move the backbone conformations!"
		close(20)
		stop
	endif

	if(recalcu_switch==0) then
		original_group=group  ! original group used for calculating RMSD
		if(flag1==1) then
			open(20, file="warning.txt", access="append")
				write(20,*) "There are Pro residues on the original peptide chain!"
				write(20,*) "please prepare a supplemental <backup4backbone.txt> file"
			close(20)
			open(20, file="backup4backbone.txt", status="old")
				do while(.true.)
					read(20, *, iostat=status) char, anum, atype, name, anum, x, y, z
					if(status.ne.0) exit
					groupdata_backup(anum)%coo(1)=x
					groupdata_backup(anum)%coo(2)=y
					groupdata_backup(anum)%coo(3)=z
				enddo
			close(20)
		endif
	else
		original_group%cnum1=0
		original_group%cnum2=0
		original_group%cnum3=0

		open(10, file=filename2)
		do while(.true.)
			read(10, *, iostat=status) char, anum, atype, name, anum, x, y, z
			if(status.ne.0) exit
			if(anum.le.sitenum) then
				if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
				   .or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
	  		  	   .or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
				   .or.name=="TYR".or.name=="VAL") then
					if (anum.eq.1) then
						name = "N" // name
					elseif (anum.eq.sitenum) then
						name = "C" // name
					endif
					original_group(anum)%gtype=name
					if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
						.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
						original_group(anum)%cnum1=original_group(anum)%cnum1+1
						original_group(anum)%atype1(original_group(anum)%cnum1)=atype
						original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
						original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
						original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
					elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
						original_group(anum)%cnum3=original_group(anum)%cnum3+1
						original_group(anum)%atype3(original_group(anum)%cnum3)=atype
						original_group(anum)%coo3(original_group(anum)%cnum3,1)=x
						original_group(anum)%coo3(original_group(anum)%cnum3,2)=y
						original_group(anum)%coo3(original_group(anum)%cnum3,3)=z
					else
						original_group(anum)%cnum2=original_group(anum)%cnum2+1
						original_group(anum)%atype2(original_group(anum)%cnum2)=atype
						original_group(anum)%coo2(original_group(anum)%cnum2,1)=x
						original_group(anum)%coo2(original_group(anum)%cnum2,2)=y
						original_group(anum)%coo2(original_group(anum)%cnum2,3)=z
					endif
				else
					open(20, file="error.txt", access="append")
						write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
						write(20,*) "Please check whether the group type is right or not!"
					close(20)
					stop
				endif

			else
				if(receptor_name=="rna") then
					if(name=="RC".or.name=="RA".or.name=="RU".or.name=="RG".or.name=="STA".or.name=="SMU".or.name=="1MA".or.name=="PSU"  &
				 		.or.name=="5MU".or.name=="RC5".or.name=="RA5".or.name=="RU5".or.name=="RG5".or.name=="RC3".or.name=="RA3" &
				   		.or.name=="RU3".or.name=="RG3".or.name=="2SU".or.name=="6TA") then
						original_group(anum)%gtype=name
						if(atype=="HO5'".or.atype=="H5T".or.atype=="P".or.atype=="O1P".or.atype=="OP1".or.atype=="O2P".or. & 
						atype=="OP2".or.atype=="O5'".or.atype=="C5'".or.atype=="1H5'".or.atype=="H5'".or.atype=="2H5'".or. & 
						atype=="H5''".or.atype=="C4'".or.atype=="H4'".or.atype=="O4'".or.atype=="C1'".or.atype=="H1'".or.atype=="H5'1" &
						.or.atype=="H5'2") then
							if(atype=="H5'1") then
								original_group(anum)%cnum1=original_group(anum)%cnum1+1
								original_group(anum)%atype1(original_group(anum)%cnum1)="1H5'"
								original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
								original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
								original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
							elseif(atype=="H5'2") then
								original_group(anum)%cnum1=original_group(anum)%cnum1+1
								original_group(anum)%atype1(original_group(anum)%cnum1)="2H5'"
								original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
								original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
								original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
							else
								original_group(anum)%cnum1=original_group(anum)%cnum1+1
								original_group(anum)%atype1(original_group(anum)%cnum1)=atype
								original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
								original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
								original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
							endif
						elseif(atype=="C3'".or.atype=="H3'".or.atype=="C2'".or.atype=="1H2'".or.atype=="H2'".or.atype=="O2'" &
							.or.atype=="2HO'".or.atype=="HO2'".or.atype=="H3T".or.atype=="HO3'".or.atype=="O3'".or.atype=="H2'1".or.atype=="HO'2") then
							if(atype=="H2'1") then
								original_group(anum)%cnum3=original_group(anum)%cnum3+1
								original_group(anum)%atype3(original_group(anum)%cnum3)="1H2'"
								original_group(anum)%coo3(original_group(anum)%cnum3,1)=x
								original_group(anum)%coo3(original_group(anum)%cnum3,2)=y
								original_group(anum)%coo3(original_group(anum)%cnum3,3)=z
							elseif(atype=="HO'2") then
								original_group(anum)%cnum3=original_group(anum)%cnum3+1
								original_group(anum)%atype3(original_group(anum)%cnum3)="2HO'"
								original_group(anum)%coo3(original_group(anum)%cnum3,1)=x
								original_group(anum)%coo3(original_group(anum)%cnum3,2)=y
								original_group(anum)%coo3(original_group(anum)%cnum3,3)=z
							else
								original_group(anum)%cnum3=original_group(anum)%cnum3+1
								original_group(anum)%atype3(original_group(anum)%cnum3)=atype
								original_group(anum)%coo3(original_group(anum)%cnum3,1)=x
								original_group(anum)%coo3(original_group(anum)%cnum3,2)=y
								original_group(anum)%coo3(original_group(anum)%cnum3,3)=z
							endif
						else
							original_group(anum)%cnum2=original_group(anum)%cnum2+1
							original_group(anum)%atype2(original_group(anum)%cnum2)=atype
							original_group(anum)%coo2(original_group(anum)%cnum2,1)=x
							original_group(anum)%coo2(original_group(anum)%cnum2,2)=y
							original_group(anum)%coo2(original_group(anum)%cnum2,3)=z
						endif
					else
						open(20, file="error.txt", access="append")
							write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
							write(20,*) "Please check whether the group type is right or not!"
						close(20)
						stop
					endif

				elseif(receptor_name=="peptide") then
					if(name=="ALA".or.name=="ARG".or.name=="ASN".or.name=="ASP".or.name=="CYS".or.name=="GLN" &
						.or.name=="GLU".or.name=="GLY".or.name=="HIE".or.name=="ILE".or.name=="LEU".or.name=="LYS" &
						.or.name=="MET".or.name=="PHE".or.name=="PRO".or.name=="SER".or.name=="THR".or.name=="TRP" &
						.or.name=="TYR".or.name=="VAL".or.name=="NALA".or.name=="NARG".or.name=="NASN".or.name=="NASP" &
						.or.name=="NCYS".or.name=="NGLN".or.name=="NGLU".or.name=="NGLY".or.name=="NHIE".or.name=="NILE" &
						.or.name=="NLEU".or.name=="NLYS".or.name=="NMET".or.name=="NPHE".or.name=="NPRO".or.name=="NSER" &
						.or.name=="NTHR".or.name=="NTRP".or.name=="NTYR".or.name=="NVAL".or.name=="CALA".or.name=="CARG" &
						.or.name=="CASN".or.name=="CASP".or.name=="CCYS".or.name=="CGLN".or.name=="CGLU".or.name=="CGLY" &
						.or.name=="CHIE".or.name=="CILE".or.name=="CLEU".or.name=="CLYS".or.name=="CMET".or.name=="CPHE" &
						.or.name=="CPRO".or.name=="CSER".or.name=="CTHR".or.name=="CTRP".or.name=="CTYR".or.name=="CVAL" &
			            .or.name=="TYX".or.name=="ARN".or.name=="HIP".or.name=="CYT".or.name=="LYN".or.name=="GLH" &
						.or.name=="ASH".or.name=="NTYX".or.name=="NARN".or.name=="NHIP".or.name=="NCYT".or.name=="NLYN" &
						.or.name=="NGLH".or.name=="NASH".or.name=="CTYX".or.name=="CARN".or.name=="CHIP".or.name=="CCYT" &
						.or.name=="CLYN".or.name=="CGLH".or.name=="CASH") then
						original_group(anum)%gtype=name
						if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
							.or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
							original_group(anum)%cnum1=original_group(anum)%cnum1+1
							original_group(anum)%atype1(original_group(anum)%cnum1)=atype
							original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
							original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
							original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
						elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
							original_group(anum)%cnum3=original_group(anum)%cnum3+1
							original_group(anum)%atype3(original_group(anum)%cnum3)=atype
							original_group(anum)%coo3(original_group(anum)%cnum3,1)=x
							original_group(anum)%coo3(original_group(anum)%cnum3,2)=y
							original_group(anum)%coo3(original_group(anum)%cnum3,3)=z
						else
							original_group(anum)%cnum2=original_group(anum)%cnum2+1
							original_group(anum)%atype2(original_group(anum)%cnum2)=atype
							original_group(anum)%coo2(original_group(anum)%cnum2,1)=x
							original_group(anum)%coo2(original_group(anum)%cnum2,2)=y
							original_group(anum)%coo2(original_group(anum)%cnum2,3)=z
						endif
					else
						open(20, file="error.txt", access="append")
							write(20,*) name, "in the file of", filename1, "is unknown in the LIB!"
							write(20,*) "Please check whether the group type is right or not!"
						close(20)
						stop
					endif

				elseif(receptor_name=="molecule") then
						original_group(anum)%gtype=name
						original_group(anum)%cnum1=original_group(anum)%cnum1+1
						original_group(anum)%atype1(original_group(anum)%cnum1)=atype
						original_group(anum)%coo1(original_group(anum)%cnum1,1)=x
						original_group(anum)%coo1(original_group(anum)%cnum1,2)=y
						original_group(anum)%coo1(original_group(anum)%cnum1,3)=z
				endif
			endif
		enddo
		close(10)
		open(5, file="backup4backbone.txt", status="old")
			do i=1, sitenum
				read(5,*) groupdata_backup(i)%coo(1), groupdata_backup(i)%coo(2), groupdata_backup(i)%coo(3)
			enddo
			read(5,*) flag4conformer
			read(5,*) energy_min
		close(5)
	endif

	return
	end subroutine readpdb

	subroutine PH_checking(group, gnum)
	implicit none
	integer*8							:: i, ic, count, gnum
	type(groupdetails)				:: group(gnum)

	count=0
	if(ph_value.le.3.9) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY".or. &
				group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" & 
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIP".or.group(ic)%gtype=="NHIP".or.group(ic)%gtype=="CHIP" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLH".or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="CGLH" &
				.or.group(ic)%gtype=="ASH".or.group(ic)%gtype=="NASH".or.group(ic)%gtype=="CASH") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY" & 
			.or.group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" & 
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR" & 
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIP".or.group(ic)%gtype=="NHIP".or.group(ic)%gtype=="CHIP" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLH".or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="CGLH" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY" &
			.or.group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET &
				".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIP".or.group(ic)%gtype=="NHIP".or.group(ic)%gtype=="CHIP" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY" & 
			.or.group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY" & 
			.or.group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY".or. &
			group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYX".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CTYX" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY".or. &
			group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYX".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CTYX" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG" &
				.or.group(ic)%gtype=="LYN".or.group(ic)%gtype=="NLYN".or.group(ic)%gtype=="CLYN" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	elseif(ph_value.ge.12.5) then
		if(receptor_name=="rna") then
			count=sitenum
		elseif(receptor_name=="peptide") then
			count=gnum
		elseif(receptor_name=="molecule") then
			count=sitenum
		endif
		do i=1, count
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY".or. &
			group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL" &
				.or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE" &
				.or.group(ic)%gtype=="TYX".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CTYX" &
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP" &
				.or.group(ic)%gtype=="ARN".or.group(ic)%gtype=="NARN".or.group(ic)%gtype=="CARN" &
				.or.group(ic)%gtype=="LYN".or.group(ic)%gtype=="NLYN".or.group(ic)%gtype=="CLYN" &
				.or.group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER" &
				.or.group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR" &
				.or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN" &
				.or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE" &
				.or.group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO" &
				.or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT" &
				.or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA" &
				.or.group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU" &
				.or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			else
				open(20, file="error.txt", access="append")
					write(20,*) group(ic)%gtype, "of the", ic, "residue could not be found in the PH of", ph_value
					write(20,*) "Please check whether the group name is right or not!"
				close(20)
				stop
			endif
		enddo

	endif

	return
	end subroutine PH_checking

	subroutine generatepdb(cycle, step, attempt, group)
	implicit none
	integer*8							:: cycle, step, attempt
	integer*8							:: i, j, atomnum
	character*5						:: stepchar, attemptchar, cyclechar
	character*4 					:: res_type
	type(groupdetails)				:: group(gnum)

	write(cyclechar,"(i3)") cycle
	write(stepchar,"(i5)") step
	write(attemptchar,"(i4)") attempt

	open(10,file='pdbfiles/'//trim(adjustl(cyclechar))//'_'//trim(adjustl(stepchar)) &
	//'_'//trim(adjustl(attemptchar))//'.pdb') 
		atomnum=1
		do i=1, gnum

			!separate peptide from receptor
			if (i.eq.sitenum+1)then
				write(10,*) "TER"
			endif

			!trim residue name for N- and C-termini
			if (i.eq.1.or.i.eq.sitenum) then
				res_type = group(i)%gtype(2:)
			else
				res_type = group(i)%gtype
			endif

			do j=1, group(i)%cnum1
				if(len_trim(group(i)%atype1(j))==4) then
					write(10,2) "ATOM", atomnum, group(i)%atype1(j), " ", res_type, i, &
					group(i)%coo1(j,1), group(i)%coo1(j,2), group(i)%coo1(j,3)
				else
					write(10,1) "ATOM", atomnum, group(i)%atype1(j), res_type, i, &
					group(i)%coo1(j,1), group(i)%coo1(j,2), group(i)%coo1(j,3)
				endif
				atomnum=atomnum+1
			enddo
			do j=1, group(i)%cnum2
				if(len_trim(group(i)%atype2(j))==4) then
					write(10,2) "ATOM", atomnum, group(i)%atype2(j), " ", res_type, i, &
					group(i)%coo2(j,1), group(i)%coo2(j,2), group(i)%coo2(j,3)
				else
					write(10,1) "ATOM", atomnum, group(i)%atype2(j), res_type, i, &
					group(i)%coo2(j,1), group(i)%coo2(j,2), group(i)%coo2(j,3)
				endif
				atomnum=atomnum+1
			enddo
			do j=1, group(i)%cnum3
				if(len_trim(group(i)%atype3(j))==4) then
					write(10,2) "ATOM", atomnum, group(i)%atype3(j), " ", res_type, i, &
					group(i)%coo3(j,1), group(i)%coo3(j,2), group(i)%coo3(j,3)
				else
					write(10,1) "ATOM", atomnum, group(i)%atype3(j), res_type, i, &
					group(i)%coo3(j,1), group(i)%coo3(j,2), group(i)%coo3(j,3)
				endif
				atomnum=atomnum+1
			enddo
1			format(a4,i7,a6,a4,i5,f12.3,2f8.3,2f6.2,a16)
2			format(a4,i7,a5,a1,a4,i5,f12.3,2f8.3,2f6.2,a16)
		enddo
	close(10)

	return
	end subroutine generatepdb

end module pdbfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mathfunction

	use constant

	contains
	subroutine normalvector(rsta, rmid, rend, r_nor)
	implicit none
	double precision						:: rsta(3), rmid(3), rend(3), r_nor(3)
	double precision						:: a(3), b(3)

	a=rsta-rmid
	b=rend-rmid

	r_nor(1)=a(2)*b(3)-a(3)*b(2)
	r_nor(2)=a(3)*b(1)-a(1)*b(3)
	r_nor(3)=a(1)*b(2)-a(2)*b(1)

	return
	end subroutine normalvector

	subroutine vectorrotation(rsta, rend, m)
	implicit none
	double precision						:: rsta(3), rend(3), m(3,3)
	double precision						:: r_cropro(3), a(3), a1(3,3), a2(3,3), a3(3,3)
	double precision						:: absrsta, absrend, r_dotpro, cos, sin

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

	r_dotpro=dot_product(rsta, rend)

	r_cropro(1)=rsta(2)*rend(3)-rsta(3)*rend(2)
	r_cropro(2)=rsta(3)*rend(1)-rsta(1)*rend(3)
	r_cropro(3)=rsta(1)*rend(2)-rsta(2)*rend(1)

	cos=r_dotpro/(absrsta*absrend)
	sin=sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))/(absrsta*absrend)

	a(1)=r_cropro(1)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(2)=r_cropro(2)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(3)=r_cropro(3)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return
	end subroutine vectorrotation

	subroutine axisrotation(a, cos, sin, m)
	implicit none
	double precision					:: cos, sin
	double precision					:: a(3), a1(3,3), a2(3,3), a3(3,3)
	double precision					:: m(3,3)

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return
	end subroutine axisrotation

	subroutine transformatrix(bond_angle, dihedral_angle, T)
	implicit none
	double precision					:: bond_angle, dihedral_angle
	double precision					:: T(3,3)

	T(1,1)=cosd(bond_angle)
	T(1,2)=sind(bond_angle)
	T(1,3)=0.0
	T(2,1)=-sind(bond_angle)*cosd(dihedral_angle)
	T(2,2)=cosd(bond_angle)*cosd(dihedral_angle)
	T(2,3)=sind(dihedral_angle)
	T(3,1)=sind(bond_angle)*sind(dihedral_angle)
	T(3,2)=-cosd(bond_angle)*sind(dihedral_angle)
	T(3,3)=cosd(dihedral_angle)

	return
	end subroutine transformatrix

	subroutine phipsiomg_angle(p1, p2, p3, p4, angle)
	implicit none
	double precision						:: p1(3), p2(3), p3(3), p4(3)
	double precision						:: angle, angle_T1, angle_T2
	double precision						:: rsta(3), rend(3), rmid(3)
	double precision						:: r_1(3), r_2(3)
	double precision						:: absrsta, absrend, r_dotpro, ratio

	call normalvector(p1, p2, p3, rend)
	call normalvector(p2, p3, p4, rsta)

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))
	r_dotpro=dot_product(rsta, rend)

	!added this logic since we could occasionally get the ratio to be greater than 1 due to floating point error, which would cause program to crash
	ratio = r_dotpro/(absrsta*absrend)
	if (abs(ratio).gt.1) then
		if (ratio.lt.0) then
			ratio = -1
		else
			ratio = 1
		endif
	endif
	angle_T1=acosd(ratio)

	if(abs(180.0-angle_T1).le.(0.5)) then
		angle=180.0
	elseif(abs(angle_T1).le.(0.4)) then
		angle=0.0
	else
		rmid=0.0
		r_2=p3-p2
		call normalvector(rsta, rmid, rend, r_1)

		absrsta=sqrt(r_1(1)*r_1(1)+r_1(2)*r_1(2)+r_1(3)*r_1(3))
		absrend=sqrt(r_2(1)*r_2(1)+r_2(2)*r_2(2)+r_2(3)*r_2(3))

		r_dotpro=dot_product(r_1, r_2)
		if(abs(r_dotpro/(absrsta*absrend)-1.0).le.(0.1)) then
			angle_T2=0.00
		elseif(abs(r_dotpro/(absrsta*absrend)+1.0).le.(0.1)) then
			angle_T2=180
		else
			angle_T2=acosd(r_dotpro/(absrsta*absrend))
		endif
		if(angle_T2.gt.90) then
			angle=-angle_T1
		else
			angle=angle_T1
		endif
	endif

	return
	end subroutine phipsiomg_angle

	subroutine rmsd_calculation(group, rmsd)
	implicit none
	integer*8							:: i, j
	double precision							:: rmsd
	double precision							:: N1(3), N2(3), CA1(3), CA2(3), C1(3), C2(3)
	type(groupdetails)				:: group(gnum)

	rmsd=0.0
	do i=1, sitenum
		do j=1, group(i)%cnum1
			if(group(i)%atype1(j)=="N") then
				N1(1)=group(i)%coo1(j,1)
				N1(2)=group(i)%coo1(j,2)
				N1(3)=group(i)%coo1(j,3)
			elseif(group(i)%atype1(j)=="CA") then
				CA1(1)=group(i)%coo1(j,1)
				CA1(2)=group(i)%coo1(j,2)
				CA1(3)=group(i)%coo1(j,3)
			endif
		enddo
		do j=1, group(i)%cnum3
			if(group(i)%atype3(j)=="C") then
				C1(1)=group(i)%coo3(j,1)
				C1(2)=group(i)%coo3(j,2)
				C1(3)=group(i)%coo3(j,3)
			endif
		enddo

		do j=1, original_group(i)%cnum1
			if(original_group(i)%atype1(j)=="N") then
				N2(1)=original_group(i)%coo1(j,1)
				N2(2)=original_group(i)%coo1(j,2)
				N2(3)=original_group(i)%coo1(j,3)
			elseif(original_group(i)%atype1(j)=="CA") then
				CA2(1)=original_group(i)%coo1(j,1)
				CA2(2)=original_group(i)%coo1(j,2)
				CA2(3)=original_group(i)%coo1(j,3)
			endif
		enddo
		do j=1, original_group(i)%cnum3
			if(original_group(i)%atype3(j)=="C") then
				C2(1)=original_group(i)%coo3(j,1)
				C2(2)=original_group(i)%coo3(j,2)
				C2(3)=original_group(i)%coo3(j,3)
			endif
		enddo
		rmsd=rmsd+(N1(1)-N2(1))**2+(N1(2)-N2(2))**2+(N1(3)-N2(3))**2+(CA1(1)-CA2(1))**2+  &
			 (CA1(2)-CA2(2))**2+(CA1(3)-CA2(3))**2+(C1(1)-C2(1))**2+(C1(2)-C2(2))**2+(C1(3)-C2(3))**2
	enddo

	rmsd=sqrt(rmsd/(3*sitenum))

	return
	end subroutine rmsd_calculation

	subroutine backbonemove_criterion(rmsd, feedback)
	implicit none
	integer*8							:: feedback
	double precision							:: rmsd

	feedback=0
	if(rmsd.le.rmsd_max.and.rmsd.gt.rmsd_min) feedback=1

	return
	end subroutine backbonemove_criterion

end module mathfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module database

	use constant
	use randomgenerator
	use mathfunction

	contains
	subroutine rotamerlib
	implicit none
	integer*8							:: status, grade, rotanum, anum, num, i, j
	double precision							:: x, y, z
	integer*8            			    :: cur_dihedrals(5)
	character*4						:: char, atype, name

	if(ph_value.le.3.9) then
		open(10, file="lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYR") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIP") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYS") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLH") then
					i=19
				elseif(name=="ASH") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYR") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIP") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYS") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLH") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYR") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIP") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYS") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYR") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIE") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYS") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYR") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIE") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYT") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYX") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYS") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIE") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYT") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYX") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARG") then
					i=9
				elseif(name=="LYN") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIE") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYT") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	elseif(ph_value.ge.12.5) then
		open(10, file="~/PepBD/lib/rotamer", status="old")
			do while(.true.)
				read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum
				if(status.ne.0) exit
				i=0
				if(name=="GLY") then
					i=1
				elseif(name=="LEU") then
					i=2
				elseif(name=="VAL") then
					i=3
				elseif(name=="ILE") then
					i=4
				elseif(name=="MET") then
					i=5
				elseif(name=="PHE") then
					i=6
				elseif(name=="TYX") then
					i=7
				elseif(name=="TRP") then
					i=8
				elseif(name=="ARN") then
					i=9
				elseif(name=="LYN") then
					i=10
				elseif(name=="SER") then
					i=11
				elseif(name=="THR") then
					i=12
				elseif(name=="ASN") then
					i=13
				elseif(name=="GLN") then
					i=14
				elseif(name=="HIE") then
					i=15
				elseif(name=="PRO") then
					i=16
				elseif(name=="CYT") then
					i=17
				elseif(name=="ALA") then
					i=18
				elseif(name=="GLU") then
					i=19
				elseif(name=="ASP") then
					i=20
				endif

				if(i.ne.0) then
					aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
					if(grade.ne.0) then
						do j=1, rotanum
							read(10,*) cur_dihedrals(1:grade)
							aa_lib(i)%dihedralangle(j,1:grade) = cur_dihedrals(1:grade)
						enddo
					endif
				endif
			enddo
		close(10)
	endif

	aa_lib%cnum1=0
	aa_lib%cnum2=0
	aa_lib%cnum3=0

	do i=1, 20
		open (10, file='~/PepBD/lib/RotamerLibrary/'//trim(aa_lib(i)%gtype), status="old")
		do while(.true.)
			read(10, *, iostat=status) char, anum, atype, name, char, num, x, y, z
			if(status.ne.0) exit
			if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
			   .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then
				aa_lib(i)%cnum1=aa_lib(i)%cnum1+1
				aa_lib(i)%atype1(aa_lib(i)%cnum1)=atype
				aa_lib(i)%coo1(aa_lib(i)%cnum1,1)=x
				aa_lib(i)%coo1(aa_lib(i)%cnum1,2)=y
				aa_lib(i)%coo1(aa_lib(i)%cnum1,3)=z
			elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
				aa_lib(i)%cnum3=aa_lib(i)%cnum3+1
				aa_lib(i)%atype3(aa_lib(i)%cnum3)=atype
				aa_lib(i)%coo3(aa_lib(i)%cnum3,1)=x
				aa_lib(i)%coo3(aa_lib(i)%cnum3,2)=y
				aa_lib(i)%coo3(aa_lib(i)%cnum3,3)=z
			else
				aa_lib(i)%cnum2=aa_lib(i)%cnum2+1
				aa_lib(i)%atype2(aa_lib(i)%cnum2)=atype
				aa_lib(i)%coo2(aa_lib(i)%cnum2,1)=x
				aa_lib(i)%coo2(aa_lib(i)%cnum2,2)=y
				aa_lib(i)%coo2(aa_lib(i)%cnum2,3)=z
			endif
		end do
		close(10)
	enddo

	return
	end subroutine rotamerlib

	subroutine findrotamer(ran_resi, group, name_original, rotanum, aa_group)
	implicit none
	integer*8								:: ran_resi, rotanum, i, j, k, l, ip
	integer*8								:: grade, grade_num(6), monitor(6)
	double precision								:: nr(3), car(3), cr(3), r_norpep(3)
	double precision								:: aa_nr(3), aa_car(3), aa_cr(3), r_norrot(3)
	double precision								:: r_nca(3), aa_r_nca(3), r_trans(3)
	double precision								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	double precision								:: delta_chi, cos_angle, sin_angle
	double precision								:: temp1(20,3), temp2(60,3), temp3(20,3)
	character*4							:: name_original

	type(groupdetails)					:: group(gnum), aa_group(40)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6)


	if(name_original=="GLY".or.name_original=="NGLY".or.name_original=="CGLY") then
		ip=1
	elseif(name_original=="LEU".or.name_original=="NLEU".or.name_original=="CLEU") then
		ip=2
	elseif(name_original=="VAL".or.name_original=="NVAL".or.name_original=="CVAL") then
		ip=3
	elseif(name_original=="ILE".or.name_original=="NILE".or.name_original=="CILE") then
		ip=4
	elseif(name_original=="MET".or.name_original=="NMET".or.name_original=="CMET") then
		ip=5
	elseif(name_original=="PHE".or.name_original=="NPHE".or.name_original=="CPHE") then
		ip=6
	elseif(name_original=="TYR".or.name_original=="NTYR".or.name_original=="CTYR".or. &
	       name_original=="TYX".or.name_original=="NTYX".or.name_original=="CTYX") then
		ip=7
	elseif(name_original=="TRP".or.name_original=="NTRP".or.name_original=="CTRP") then
		ip=8
	elseif(name_original=="ARG".or.name_original=="NARG".or.name_original=="CARG".or. &
	       name_original=="ARN".or.name_original=="NARN".or.name_original=="CARN") then
		ip=9
	elseif(name_original=="LYS".or.name_original=="NLYS".or.name_original=="CLYS".or. &
	       name_original=="LYN".or.name_original=="NLYN".or.name_original=="CLYN") then
		ip=10
	elseif(name_original=="SER".or.name_original=="NSER".or.name_original=="CSER") then
		ip=11
	elseif(name_original=="THR".or.name_original=="NTHR".or.name_original=="CTHR") then
		ip=12
	elseif(name_original=="ASN".or.name_original=="NASN".or.name_original=="CASN") then
		ip=13
	elseif(name_original=="GLN".or.name_original=="NGLN".or.name_original=="CGLN") then
		ip=14
	elseif(name_original=="HIE".or.name_original=="NHIE".or.name_original=="CHIE".or. &
	       name_original=="HIP".or.name_original=="NHIP".or.name_original=="CHIP") then
		ip=15
	elseif(name_original=="PRO".or.name_original=="NPRO".or.name_original=="CPRO") then
		ip=16
	elseif(name_original=="CYS".or.name_original=="NCYS".or.name_original=="CCYS".or. &
	       name_original=="CYT".or.name_original=="NCYT".or.name_original=="CCYT") then
		ip=17
	elseif(name_original=="ALA".or.name_original=="NALA".or.name_original=="CALA") then
		ip=18
	elseif(name_original=="GLU".or.name_original=="NGLU".or.name_original=="CGLU".or. &
	       name_original=="GLH".or.name_original=="NGLH".or.name_original=="CGLH") then
		ip=19
	elseif(name_original=="ASP".or.name_original=="NASP".or.name_original=="CASP".or. &
	       name_original=="ASH".or.name_original=="NASH".or.name_original=="CASH") then
		ip=20
	endif

	rotanum=aa_lib(ip)%rotanum
	aa_group(1)%cnum1=aa_lib(ip)%cnum1; aa_group(1)%cnum2=aa_lib(ip)%cnum2; aa_group(1)%cnum3=aa_lib(ip)%cnum3
	aa_group(1)%gtype=name_original
	aa_group(1)%atype1=aa_lib(ip)%atype1; aa_group(1)%atype2=aa_lib(ip)%atype2; aa_group(1)%atype3=aa_lib(ip)%atype3
	aa_group(1)%coo1=aa_lib(ip)%coo1; aa_group(1)%coo2=aa_lib(ip)%coo2; aa_group(1)%coo3=aa_lib(ip)%coo3

	do i=1, group(ran_resi)%cnum1
		if(group(ran_resi)%atype1(i)=="N") then
			nr(1)=group(ran_resi)%coo1(i,1)
			nr(2)=group(ran_resi)%coo1(i,2)
			nr(3)=group(ran_resi)%coo1(i,3)
		elseif(group(ran_resi)%atype1(i)=="CA") then
			car(1)=group(ran_resi)%coo1(i,1)
			car(2)=group(ran_resi)%coo1(i,2)
			car(3)=group(ran_resi)%coo1(i,3)
		endif
	enddo
	do i=1, group(ran_resi)%cnum3
		if(group(ran_resi)%atype3(i)=="C") then
			cr(1)=group(ran_resi)%coo3(i,1)
			cr(2)=group(ran_resi)%coo3(i,2)
			cr(3)=group(ran_resi)%coo3(i,3)
		endif
	enddo

	call normalvector(nr, car, cr, r_norpep)

	r_nca(1)=nr(1)-car(1)
	r_nca(2)=nr(2)-car(2)
	r_nca(3)=nr(3)-car(3)

	do i=1, aa_group(1)%cnum1
		if(aa_group(1)%atype1(i)=="N") then
			aa_nr(1)=aa_group(1)%coo1(i,1)
			aa_nr(2)=aa_group(1)%coo1(i,2)
			aa_nr(3)=aa_group(1)%coo1(i,3)
		elseif(aa_group(1)%atype1(i)=="CA") then
			aa_car(1)=aa_group(1)%coo1(i,1)
			aa_car(2)=aa_group(1)%coo1(i,2)
			aa_car(3)=aa_group(1)%coo1(i,3)
		endif
	enddo
	do i=1, aa_group(1)%cnum3
		if(aa_group(1)%atype3(i)=="C") then
			aa_cr(1)=aa_group(1)%coo3(i,1)
			aa_cr(2)=aa_group(1)%coo3(i,2)
			aa_cr(3)=aa_group(1)%coo3(i,3)
		endif
	enddo

	call normalvector(aa_nr, aa_car, aa_cr, r_norrot)

	call vectorrotation(r_norrot, r_norpep, m)

	temp1=matmul(aa_group(1)%coo1, m)
	aa_group(1)%coo1=temp1

	temp2=matmul(aa_group(1)%coo2, m)
	aa_group(1)%coo2=temp2

	temp3=matmul(aa_group(1)%coo3, m)
	aa_group(1)%coo3=temp3

	do i=1, aa_group(1)%cnum1
		if(aa_group(1)%atype1(i)=="N") then
			aa_nr(1)=aa_group(1)%coo1(i,1)
			aa_nr(2)=aa_group(1)%coo1(i,2)
			aa_nr(3)=aa_group(1)%coo1(i,3)
		elseif(aa_group(1)%atype1(i)=="CA") then
			aa_car(1)=aa_group(1)%coo1(i,1)
			aa_car(2)=aa_group(1)%coo1(i,2)
			aa_car(3)=aa_group(1)%coo1(i,3)
		endif
	enddo

	aa_r_nca(1)=aa_nr(1)-aa_car(1)
	aa_r_nca(2)=aa_nr(2)-aa_car(2)
	aa_r_nca(3)=aa_nr(3)-aa_car(3)

	call vectorrotation(aa_r_nca, r_nca, m)

	temp1=matmul(aa_group(1)%coo1, m)
	aa_group(1)%coo1=temp1

	temp2=matmul(aa_group(1)%coo2, m)
	aa_group(1)%coo2=temp2

	temp3=matmul(aa_group(1)%coo3, m)
	aa_group(1)%coo3=temp3

	do i=1, aa_group(1)%cnum1
		if(aa_group(1)%atype1(i)=="CA") then
			aa_car(1)=aa_group(1)%coo1(i,1)
			aa_car(2)=aa_group(1)%coo1(i,2)
			aa_car(3)=aa_group(1)%coo1(i,3)
		endif
	enddo

	r_trans(1)=car(1)-aa_car(1)
	r_trans(2)=car(2)-aa_car(2)
	r_trans(3)=car(3)-aa_car(3)

	do i=1, aa_group(1)%cnum1
		aa_group(1)%coo1(i,1)=anint((aa_group(1)%coo1(i,1)+r_trans(1))*1000)/1000
		aa_group(1)%coo1(i,2)=anint((aa_group(1)%coo1(i,2)+r_trans(2))*1000)/1000
		aa_group(1)%coo1(i,3)=anint((aa_group(1)%coo1(i,3)+r_trans(3))*1000)/1000
	enddo
	do i=1, aa_group(1)%cnum2
		aa_group(1)%coo2(i,1)=anint((aa_group(1)%coo2(i,1)+r_trans(1))*1000)/1000
		aa_group(1)%coo2(i,2)=anint((aa_group(1)%coo2(i,2)+r_trans(2))*1000)/1000
		aa_group(1)%coo2(i,3)=anint((aa_group(1)%coo2(i,3)+r_trans(3))*1000)/1000
	enddo
	do i=1, aa_group(1)%cnum3
		aa_group(1)%coo3(i,1)=anint((aa_group(1)%coo3(i,1)+r_trans(1))*1000)/1000
		aa_group(1)%coo3(i,2)=anint((aa_group(1)%coo3(i,2)+r_trans(2))*1000)/1000
		aa_group(1)%coo3(i,3)=anint((aa_group(1)%coo3(i,3)+r_trans(3))*1000)/1000
	enddo

	if(rotanum.le.1) goto 10

	grade_num=0
	if(aa_group(1)%gtype=="VAL".or.aa_group(1)%gtype=="NVAL".or.aa_group(1)%gtype=="CVAL") then
		grade=1
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1)
				Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2)
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo
	elseif(aa_group(1)%gtype=="LEU".or.aa_group(1)%gtype=="NLEU".or.aa_group(1)%gtype=="CLEU") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="ILE".or.aa_group(1)%gtype=="NILE".or.aa_group(1)%gtype=="CILE") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB".or.aa_group(1)%atype2(i)=="CG2".or.aa_group(1)%atype2(i)=="HG21".or. &
				aa_group(1)%atype2(i)=="HG22".or.aa_group(1)%atype2(i)=="HG23".or.aa_group(1)%atype2(i)=="CG1") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG1") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="PHE".or.aa_group(1)%gtype=="NPHE".or.aa_group(1)%gtype=="CPHE") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="TRP".or.aa_group(1)%gtype=="NTRP".or.aa_group(1)%gtype=="CTRP") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="TYR".or.aa_group(1)%gtype=="NTYR".or.aa_group(1)%gtype=="CTYR".or. &
	       aa_group(1)%gtype=="TYX".or.aa_group(1)%gtype=="NTYX".or.aa_group(1)%gtype=="CTYX") then
		grade=3
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HH") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="OH") monitor(3)=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="SER".or.aa_group(1)%gtype=="NSER".or.aa_group(1)%gtype=="CSER") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="OG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="OG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="THR".or.aa_group(1)%gtype=="NTHR".or.aa_group(1)%gtype=="CTHR") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HG1") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="OG1") monitor(2)=grade_num(2)
			endif
		enddo
	elseif(aa_group(1)%gtype=="CYS".or.aa_group(1)%gtype=="NCYS".or.aa_group(1)%gtype=="CCYS".or. &
	       aa_group(1)%gtype=="CYT".or.aa_group(1)%gtype=="NCYT".or.aa_group(1)%gtype=="CCYT") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="SG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="SG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="MET".or.aa_group(1)%gtype=="NMET".or.aa_group(1)%gtype=="CMET") then
		grade=4
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2);
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="SD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="SD") monitor(3)=grade_num(3)
			elseif(aa_group(1)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=aa_group(1)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=aa_group(1)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(aa_group(1)%gtype=="ASN".or.aa_group(1)%gtype=="NASN".or.aa_group(1)%gtype=="CASN") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="GLN".or.aa_group(1)%gtype=="NGLN".or.aa_group(1)%gtype=="CGLN") then
		grade=3
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(aa_group(1)%gtype=="ASP".or.aa_group(1)%gtype=="NASP".or.aa_group(1)%gtype=="CASP".or. &
	       aa_group(1)%gtype=="ASH".or.aa_group(1)%gtype=="NASH".or.aa_group(1)%gtype=="CASH") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="GLU".or.aa_group(1)%gtype=="NGLU".or.aa_group(1)%gtype=="CGLU".or. &
	       aa_group(1)%gtype=="GLH".or.aa_group(1)%gtype=="NGLH".or.aa_group(1)%gtype=="CGLH") then
		grade=3
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(aa_group(1)%gtype=="HIE".or.aa_group(1)%gtype=="NHIE".or.aa_group(1)%gtype=="CHIE".or. &
	       aa_group(1)%gtype=="HIP".or.aa_group(1)%gtype=="NHIP".or.aa_group(1)%gtype=="CHIP") then
		grade=2
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(aa_group(1)%gtype=="LYS".or.aa_group(1)%gtype=="NLYS".or.aa_group(1)%gtype=="CLYS".or. &
	       aa_group(1)%gtype=="LYN".or.aa_group(1)%gtype=="NLYN".or.aa_group(1)%gtype=="CLYN") then
		grade=5
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(aa_group(1)%atype2(i)=="HD2".or.aa_group(1)%atype2(i)=="HD3".or.aa_group(1)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(aa_group(1)%atype2(i)=="CE") monitor(4)=grade_num(4)
			elseif(aa_group(1)%atype2(i)=="HE2".or.aa_group(1)%atype2(i)=="HE3".or.aa_group(1)%atype2(i)=="NZ") then
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=aa_group(1)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=aa_group(1)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
				if(aa_group(1)%atype2(i)=="NZ") monitor(5)=grade_num(5)
			else
				grade_num(6)=grade_num(6)+1
				Iclass(6)%member(grade_num(6),1)=aa_group(1)%coo2(i,1); Iclass(6)%member(grade_num(6),2)=aa_group(1)%coo2(i,2); 
				Iclass(6)%member(grade_num(6),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=6; index(i)%member_No=grade_num(6)
			endif
		enddo
	elseif(aa_group(1)%gtype=="ARG".or.aa_group(1)%gtype=="NARG".or.aa_group(1)%gtype=="CARG".or. &
	       aa_group(1)%gtype=="ARN".or.aa_group(1)%gtype=="NARN".or.aa_group(1)%gtype=="CARN") then
		grade=4
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(aa_group(1)%atype2(i)=="HD2".or.aa_group(1)%atype2(i)=="HD3".or.aa_group(1)%atype2(i)=="NE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(aa_group(1)%atype2(i)=="NE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=aa_group(1)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=aa_group(1)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(aa_group(1)%gtype=="PRO".or.aa_group(1)%gtype=="NPRO".or.aa_group(1)%gtype=="CPRO") then
		grade=3
		do i=1, aa_group(1)%cnum2
			if(aa_group(1)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=aa_group(1)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=aa_group(1)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(aa_group(1)%atype2(i)=="HB2".or.aa_group(1)%atype2(i)=="HB3".or.aa_group(1)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=aa_group(1)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=aa_group(1)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(aa_group(1)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(aa_group(1)%atype2(i)=="HG2".or.aa_group(1)%atype2(i)=="HG3".or.aa_group(1)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=aa_group(1)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=aa_group(1)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(aa_group(1)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=aa_group(1)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=aa_group(1)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=aa_group(1)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	endif

	if(grade.ne.aa_lib(ip)%grade) then
		open(10, file="error.txt", access="append")
			write(10,*) aa_lib(ip)%gtype
			write(10,*) "grade=", grade
			write(10,*) "aa_lib(",ip,")%grade=", aa_lib(ip)%grade
			write(10,*) "They are not equal with each other!"
		close(10)
		stop
	endif

	do i=1, aa_group(1)%cnum1
		if(aa_group(1)%atype1(i)=="CA") then
			CA(1)=aa_group(1)%coo1(i,1); CA(2)=aa_group(1)%coo1(i,2); CA(3)=aa_group(1)%coo1(i,3)
		endif
	enddo

	do i=2, rotanum
		aa_group(i)%cnum1=aa_group(1)%cnum1; aa_group(i)%cnum2=aa_group(1)%cnum2; aa_group(i)%cnum3=aa_group(1)%cnum3
		aa_group(i)%gtype=name_original
		aa_group(i)%atype1=aa_group(1)%atype1; aa_group(i)%atype2=aa_group(1)%atype2; aa_group(i)%atype3=aa_group(1)%atype3
		aa_group(i)%coo1=aa_group(1)%coo1; aa_group(i)%coo2=aa_group(1)%coo2; aa_group(i)%coo3=aa_group(1)%coo3

		Tclass=Iclass
		do j=1, grade
			delta_chi=DBLE(aa_lib(ip)%dihedralangle(i,j)-aa_lib(ip)%dihedralangle(1,j))
			cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
			if(j==1) then
				rotaxis_x=Tclass(j)%member(monitor(j),1)-CA(1)
				rotaxis_y=Tclass(j)%member(monitor(j),2)-CA(2)
				rotaxis_z=Tclass(j)%member(monitor(j),3)-CA(3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			else
				rotaxis_x=Tclass(j)%member(monitor(j),1)-Tclass(j-1)%member(monitor(j-1),1)
				rotaxis_y=Tclass(j)%member(monitor(j),2)-Tclass(j-1)%member(monitor(j-1),2)
				rotaxis_z=Tclass(j)%member(monitor(j),3)-Tclass(j-1)%member(monitor(j-1),3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			endif

			call axisrotation(rotaxis, cos_angle, sin_angle, m)

			do l=(j+1), (grade+1)
				do k=1, grade_num(l)
					Tclass(l)%member(k,1)=Tclass(l)%member(k,1)-Tclass(j)%member(monitor(j),1)
					Tclass(l)%member(k,2)=Tclass(l)%member(k,2)-Tclass(j)%member(monitor(j),2)
					Tclass(l)%member(k,3)=Tclass(l)%member(k,3)-Tclass(j)%member(monitor(j),3)
				enddo

				Tmember=matmul(Tclass(l)%member, m)
				Tclass(l)%member=Tmember

				do k=1, grade_num(l)
					Tclass(l)%member(k,1)=anint((Tclass(l)%member(k,1)+Tclass(j)%member(monitor(j),1))*1000)/1000
					Tclass(l)%member(k,2)=anint((Tclass(l)%member(k,2)+Tclass(j)%member(monitor(j),2))*1000)/1000
					Tclass(l)%member(k,3)=anint((Tclass(l)%member(k,3)+Tclass(j)%member(monitor(j),3))*1000)/1000
				enddo
			enddo
		enddo

		do l=1, aa_group(i)%cnum2
			aa_group(i)%coo2(l,1)=Tclass(index(l)%class_No)%member(index(l)%member_No,1)
			aa_group(i)%coo2(l,2)=Tclass(index(l)%class_No)%member(index(l)%member_No,2)
			aa_group(i)%coo2(l,3)=Tclass(index(l)%class_No)%member(index(l)%member_No,3)
		enddo
	enddo
10	continue

	return
	end subroutine findrotamer

	subroutine ramachandranmap
	implicit none
	integer*8					        :: i, j

	open(10, file="~/PepBD/lib/rama_map", status="old")
		do j=-179, 180
			read(10,"(360i2)") (rama_map(i,j), i=-179, 180)
		end do
	close(10)

	flavoredregion_number=0

	do i=-179, 180
		do j=-179, 180
			if(rama_map(i,j)==1.or.rama_map(i,j)==2) then
				flavoredregion_number=flavoredregion_number+1
				flavored_region(flavoredregion_number,1)=DBLE(i)
				flavored_region(flavoredregion_number,2)=DBLE(j)
			endif
		enddo
	enddo

	return
	end subroutine ramachandranmap

	subroutine energy_parameter(group, group_para, gnum)
  !This subroutine reads in the force field parameters from in an
  !input file stored in Lib/ForceField.
	implicit none
	integer*8							:: i, j, status, atomid, gnum
	double precision							:: charge, epsion, r, rborn, fs, dielecons
	character*4						:: lbres, igraph
	type(groupdetails)				:: group(gnum)
	type(energyparameters)			:: group_para(gnum)
	type(energyparameters)       ::	Tgroup_para(gnum)

	do i=1, gnum
		open(10, file='~/PepBD/lib/ForceField/'//trim(group(i)%gtype), status="old")
		read(10, *)
		do while(.true.)
			read(10, 20, iostat=status) lbres, igraph, charge, epsion, r, rborn, fs, dielecons, atomid
			if(status.ne.0) goto 30
			do j=1, group(i)%cnum1
				if(group(i)%atype1(j)==igraph) then
					Tgroup_para(i)%charge1(j)=charge
					Tgroup_para(i)%epsion1(j)=epsion
					Tgroup_para(i)%r1(j)=r
					Tgroup_para(i)%rborn1(j)=rborn
					Tgroup_para(i)%fs1(j)=fs
					Tgroup_para(i)%dielecons1(j)=dielecons
					Tgroup_para(i)%atomid1(j)=atomid
					goto 40
				endif
			enddo
			do j=1, group(i)%cnum2
				if(group(i)%atype2(j)==igraph) then
					Tgroup_para(i)%charge2(j)=charge
					Tgroup_para(i)%epsion2(j)=epsion
					Tgroup_para(i)%r2(j)=r
					Tgroup_para(i)%rborn2(j)=rborn
					Tgroup_para(i)%fs2(j)=fs
					Tgroup_para(i)%dielecons2(j)=dielecons
					Tgroup_para(i)%atomid2(j)=atomid
					goto 40
				endif
			enddo
			do j=1, group(i)%cnum3
				if(group(i)%atype3(j)==igraph) then
					Tgroup_para(i)%charge3(j)=charge
					Tgroup_para(i)%epsion3(j)=epsion
					Tgroup_para(i)%r3(j)=r
					Tgroup_para(i)%rborn3(j)=rborn
					Tgroup_para(i)%fs3(j)=fs
					Tgroup_para(i)%dielecons3(j)=dielecons
					Tgroup_para(i)%atomid3(j)=atomid
					goto 40
				endif
			enddo
40			continue
		enddo
30		continue
		close(10)

20	format(2a4, 6e16.8, i8)
	enddo

	do i=1, gnum
		do j=1, group(i)%cnum1
			if(Tgroup_para(i)%dielecons1(j)<=0.1) then
				open(10, file="error.txt", access="append")
				write(10,*) "Atom", group(i)%atype1(j), "in", group(i)%gtype, "has bad force field parameter", Tgroup_para(i)%dielecons2(j)
				write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
				close(10)
				stop
			endif
		enddo
		do j=1, group(i)%cnum2
			if(Tgroup_para(i)%dielecons2(j)<=0.1) then
				open(10, file="error.txt", access="append")
				write(10,*) "Atom", group(i)%atype2(j), "in", group(i)%gtype, "has bad force field parameter", Tgroup_para(i)%dielecons2(j)
				write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
				close(10)
				stop
			endif
		enddo
		do j=1, group(i)%cnum3
			if(Tgroup_para(i)%dielecons3(j)<=0.1) then
				open(10, file="error.txt", access="append")
				write(10,*) "Atom", group(i)%atype2(j), "in", group(i)%gtype, "has bad force field parameter", Tgroup_para(i)%dielecons2(j)
				write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
				close(10)
				stop
			endif
		enddo
	enddo
	group_para=Tgroup_para

	return
	end subroutine energy_parameter

	subroutine atom_links(group, numex, inb, numex4, inb4)
	implicit none

	integer*8				:: i, j, status, i1, j1, i2, j2, k, atomid, natom, linkindex(4)
	integer*8				:: id, linknum
	integer*8				:: ipres(gnum), numex(atom_num), inb(atom_num,20), numex4(atom_num), inb4(atom_num, 60)
	type(groupdetails)		:: group(gnum)
	type(atomlink)			:: atom(atom_num)

	natom=0

	!first read in the atom_links info for each residue in the peptide
	do i=1, sitenum !only need to create atom_links for peptide (atom_links control whether van der Waals or Coloumbic interactions measured, and these are never calculated for receptor-receptor atom pairs)
		ipres(i)=natom
		open(10, file='~/PepBD/lib/Atomlink/'//trim(group(i)%gtype), status="old")
		do while(.true.)
			read(10, 20, iostat=status) id, linknum, (linkindex(j), j=1, linknum)
			if(status.ne.0) goto 30
			atomid=ipres(i)+id
			atom(atomid)%linknum=linknum
			do j=1, linknum
				atom(atomid)%linkindex(j)=ipres(i)+linkindex(j)
			enddo
		enddo
30		continue
		close(10)
		natom=atomid
	enddo
20  format(i6, i7, 4i3)

	do i1=1, natom
		numex(i1)=0
		do j1=1, atom(i1)%linknum
			numex(i1)=numex(i1)+1
			inb(i1,numex(i1))=atom(i1)%linkindex(j1)
		enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 40
				do k=1, atom(i1)%linknum
					if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
				enddo
				numex(i1)=numex(i1)+1
				inb(i1,numex(i1))=atom(i2)%linkindex(j2)
40				continue
			enddo
		enddo
	enddo

	do i1=1, natom
		numex4(i1)=0
		do j1=1, numex(i1)
			i2=inb(i1,j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 50
				do k=1, numex(i1)
					if(atom(i2)%linkindex(j2).eq.inb(i1,k)) goto 50
				enddo
				numex4(i1)=numex4(i1)+1
				inb4(i1,numex4(i1))=atom(i2)%linkindex(j2)
50			continue
			enddo
		enddo
	enddo

	return
	end subroutine atom_links

	subroutine atom_links4sidechain(ic, group, numex, inb, numex4, inb4)
	implicit none

	integer*8						:: j, status, i1, j1, j2, k, atomid, natom, ic
	integer*8						:: linkindex(4)
	integer*8						:: id, linknum, i2
	integer*8						:: inb(60,20), inb4(60, 60), numex(60), numex4(60)
	type(groupdetails)				:: group(gnum)
	type(atomlink)					:: atom(60)

	open(10, file='~/PepBD/lib/Atomlink/'//trim(group(ic)%gtype), status="old")
	do while(.true.)
		read(10, 20, iostat=status) id, linknum, (linkindex(j), j=1, linknum)
		if(status.ne.0) goto 30
		atomid=id
		atom(atomid)%linknum=linknum
		do j=1, linknum
			atom(atomid)%linkindex(j)=linkindex(j)
		enddo
	enddo
30	continue
	close(10)
	natom=atomid
20  format(i6, i7, 4i3)

    !This obtains atoms that directly bonded to every atom in the residue
do i1=1, natom
	numex(i1)=0
	do j1=1, atom(i1)%linknum
		numex(i1)=numex(i1)+1
		inb(i1,numex(i1))=atom(i1)%linkindex(j1)
	enddo

	!This obtains atoms that two bonds away from every atom in the residue
	do j1=1, atom(i1)%linknum
		i2=atom(i1)%linkindex(j1)
		if((i2.ge.1).and.(i2.le.natom)) then
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 40
				do k=1, atom(i1)%linknum
					if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
				enddo
				numex(i1)=numex(i1)+1
				inb(i1,numex(i1))=atom(i2)%linkindex(j2)
40				continue
			enddo
		endif
	enddo
enddo

!This obtains obtains that are three bonds away from every atom in the residue (added to inb4/numex4 now, rather than inb/numex)
do i1=1, natom
	numex4(i1)=0
	do j1=1, numex(i1)
		i2=inb(i1,j1)
		if((i2.ge.1).and.(i2.le.natom)) then
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 50
				do k=1, numex(i1)
					if(atom(i2)%linkindex(j2).eq.inb(i1,k)) goto 50
				enddo
				numex4(i1)=numex4(i1)+1
				inb4(i1,numex4(i1))=atom(i2)%linkindex(j2)
50			continue
			enddo
		endif
	enddo
enddo

	return
	end subroutine atom_links4sidechain

	subroutine mc_choose_aminoacid(ic, group, aminoacid_name, res_type_old, res_type_new, trp_delta)
	implicit none
	integer*8						:: ic, ip, ip1, i, trp_delta
	integer*8						:: res_type_old, res_type_new, res_avail_count, cutoffs(6)
	double precision						:: ran2
	character*4					:: aminoacid_name
	character*3					:: res_avail(20)
	type(groupdetails)			:: group(gnum)

	if(ph_value.le.3.9) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
		       group(ic)%gtype=="GLH".or.group(ic)%gtype=="ASH") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
		       group(ic)%gtype=="NGLH".or.group(ic)%gtype=="NASH") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
		       group(ic)%gtype=="CGLH".or.group(ic)%gtype=="CASH") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
			ip=52
		endif

10		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TYR"
			elseif(ip1.eq.7) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="ARG"
			elseif(ip1.eq.2) then
				aminoacid_name="LYS"
			elseif(ip1.eq.3) then
				aminoacid_name="HIP"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="GLH"
			elseif(ip1.eq.6) then
				aminoacid_name="ASH"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			elseif(ip1.eq.3) then
				aminoacid_name="CYS"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NARG"
			elseif(ip1.eq.2) then
				aminoacid_name="NLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="NHIP"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NGLH"
			elseif(ip1.eq.6) then
				aminoacid_name="NASH"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			elseif(ip1.eq.3) then
				aminoacid_name="NCYS"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CARG"
			elseif(ip1.eq.2) then
				aminoacid_name="CLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="CHIP"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CGLH"
			elseif(ip1.eq.6) then
				aminoacid_name="CASH"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="CALA"
			elseif(ip1.eq.3) then
				aminoacid_name="CCYS"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 10

	elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
		       group(ic)%gtype=="GLH") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="ASP") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
		       group(ic)%gtype=="NGLH") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NASP") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
		       group(ic)%gtype=="CGLH") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CASP") then
			ip=62
		endif

20		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TYR"
			elseif(ip1.eq.7) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="ARG"
			elseif(ip1.eq.2) then
				aminoacid_name="LYS"
			elseif(ip1.eq.3) then
				aminoacid_name="HIP"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="GLH"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			elseif(ip1.eq.3) then
				aminoacid_name="CYS"
			endif
		elseif(ip.eq.60) then
			aminoacid_name="ASP"
			goto 100
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NARG"
			elseif(ip1.eq.2) then
				aminoacid_name="NLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="NHIP"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NGLH"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			elseif(ip1.eq.3) then
				aminoacid_name="NCYS"
			endif
		elseif(ip.eq.61) then
			aminoacid_name="NASP"
			goto 100
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CARG"
			elseif(ip1.eq.2) then
				aminoacid_name="CLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="CHIP"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CGLH"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="CALA"
			elseif(ip1.eq.3) then
				aminoacid_name="CCYS"
			endif
		elseif(ip.eq.62) then
			aminoacid_name="CASP"
			goto 100
		endif
		if(aminoacid_name==group(ic)%gtype) goto 20

	elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP") then
			ip=62
		endif

30		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TYR"
			elseif(ip1.eq.7) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="ARG"
			elseif(ip1.eq.2) then
				aminoacid_name="LYS"
			elseif(ip1.eq.3) then
				aminoacid_name="HIP"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
!			ip1=int(ran2*3-1.0e-3)+1
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="CYS"
			endif
		elseif(ip.eq.60) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="GLU"
			elseif(ip1.eq.2) then
				aminoacid_name="ASP"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NARG"
			elseif(ip1.eq.2) then
				aminoacid_name="NLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="NHIP"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
!			ip1=int(ran2*3-1.0e-3)+1
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="NCYS"
			endif
		elseif(ip.eq.61) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="NASP"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CARG"
			elseif(ip1.eq.2) then
				aminoacid_name="CLYS"
			elseif(ip1.eq.3) then
				aminoacid_name="CHIP"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
!			ip1=int(ran2*3-1.0e-3)+1
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="CALA"
!			elseif(ip1.eq.3) then
!				aminoacid_name="CCYS"
			endif
		elseif(ip.eq.62) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="CASP"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 30

	elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then

		!determine what type of amino acid is being replaced
		if (group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
			res_type_old = 1

		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU".or. &
		group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL".or. &
		group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE".or. &
		group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET".or. &
		group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE".or. &
		group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR".or. &
		group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP") then
			res_type_old = 2

		elseif (group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG".or. &
		group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS") then
			res_type_old = 3

		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER".or. &
		group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR".or. &
		group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN".or. &
		group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN".or. &
		group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE") then
			res_type_old = 4

		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="NPRO".or.group(ic)%gtype=="CPRO".or. &
		group(ic)%gtype=="ALA".or.group(ic)%gtype=="NALA".or.group(ic)%gtype=="CALA".or. &
		group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS") then
			res_type_old = 5

		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU".or. &
		group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
			res_type_old = 6

		else
			open(20, file="error.txt", access="append")
				write(20,*) "Residue name not recognized during residue replacement, terminating program"
				close(20)
			stop
		endif
		
	!determine which amino acid types can replace the selected amino acid
		res_type_new = 0
		if (below_hyd(res_type_old).eq.0) then
			!Have to replace the amino acid with one of the same type to stay within hydration range
			if (res_type_old.eq.1) then
				res_avail(1) = "GLY"
				res_avail_count = 1
				res_type_new = 1

			elseif (res_type_old.eq.2) then
				!need an extra if statement to not go over the tryptophan limit
				if (trp_count.ge.trp_limit.and.(group(ic)%gtype.ne."NTRP".or.group(ic)%gtype.ne."TRP".or.group(ic)%gtype.ne."CTRP")) then
					res_avail(1:6) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE"/)
					res_avail_count = 6
					res_type_new = 2
				else
					res_avail(1:7) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE", "TRP"/)
					res_avail_count = 7
					res_type_new = 2
				endif

			elseif (res_type_old.eq.3) then
				res_avail(1:2) = (/"ARG", "LYS"/)
				res_avail_count = 2
				res_type_new = 3

			elseif (res_type_old.eq.4) then
				res_avail(1:5) = (/"SER", "THR", "ASN", "GLN", "HIE"/)
				res_avail_count = 5
				res_type_new = 4

			elseif (res_type_old.eq.5) then
				res_avail(1) = "ALA"
				res_avail_count = 1
				res_type_new = 5

			else 
				res_avail(1:2) = (/"ASP", "GLU"/)
				res_avail_count = 2
				res_type_new = 6
			endif
		
		else
		!see which amino acid types we can increase while staying in the hydration range
			res_avail_count = 0  
			cutoffs = 0
			if (above_hyd(1).lt.0) then
				res_avail(res_avail_count+1) = ("GLY")
				res_avail_count = res_avail_count + 1
				cutoffs(1) = res_avail_count
			endif

			if (above_hyd(2).lt.0) then
				if (trp_count.ge.trp_limit.and.(group(ic)%gtype.ne."NTRP".or.group(ic)%gtype.ne."TRP".or.group(ic)%gtype.ne."CTRP")) then
					res_avail(res_avail_count+1:res_avail_count+6) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE"/)
					res_avail_count = res_avail_count + 6
				else
					res_avail(res_avail_count+1:res_avail_count+7) = (/"LEU", "ILE", "MET", "VAL", "TYR", "PHE", "TRP"/)
					res_avail_count = res_avail_count + 7
				endif
				cutoffs(2) = res_avail_count
			endif

			if (above_hyd(3).lt.0) then
				res_avail(res_avail_count+1:res_avail_count+2) = (/"ARG", "LYS"/)
				res_avail_count = res_avail_count + 2
				cutoffs(3) = res_avail_count
			endif

			if (above_hyd(4).lt.0) then
				res_avail(res_avail_count+1:res_avail_count+5) = (/"SER", "THR", "ASN", "GLN", "HIE"/)
				res_avail_count = res_avail_count + 5
				cutoffs(4) = res_avail_count
			endif

			if (above_hyd(5).lt.0) then
				res_avail(res_avail_count+1) = "ALA"
				res_avail_count = res_avail_count + 1
				cutoffs(5) = res_avail_count
			endif

			if (above_hyd(6).lt.0) then
				res_avail(res_avail_count+1:res_avail_count+2) = (/"ASP", "GLU"/)
				res_avail_count = res_avail_count + 2
				cutoffs(6) = res_avail_count
			endif

		endif

		!now pick one of the possible amino acids randomly
		call ran_gen(ran2,zero)
		ip1=min(int(ran2*res_avail_count-1.0e-3)+1, res_avail_count)
		aminoacid_name = res_avail(ip1)
		
		!Check if tryptophan count is changing
		trp_delta = 0 
		if (group(ic)%gtype.eq."TRP".or.group(ic)%gtype.eq."NTRP".or.group(ic)%gtype.eq."CTRP") then
			if (aminoacid_name.ne."TRP") trp_delta = -1
		else 
			if (aminoacid_name.eq."TRP") trp_delta = 1
		endif

		!if we are at one of the peptide termini, then add an N or C to amino acid name
		if (ic.eq.1) then
			aminoacid_name = "N" // aminoacid_name
		elseif (ic.eq.sitenum) then
			aminoacid_name = "C" // aminoacid_name
		endif

		!figure out the type of the replacement residue, if not already known
		if (res_type_new.eq.0) then
			i = 1
			do while (res_type_new.eq.0)
				if (cutoffs(i).ge.ip1) then
					res_type_new = i
				endif
				i = i + 1 
			enddo
		endif 

	elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
			group(ic)%gtype=="HIE") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
			group(ic)%gtype=="NHIE") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="NCYT") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
			group(ic)%gtype=="CHIE") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT") then
			ip=62
		endif

50		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TYR"
			elseif(ip1.eq.7) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="ARG"
			elseif(ip1.eq.2) then
				aminoacid_name="LYS"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="HIE"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			endif
		elseif(ip.eq.60) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="GLU"
			elseif(ip1.eq.2) then
				aminoacid_name="ASP"
			elseif(ip1.eq.3) then
				aminoacid_name="CYT"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NARG"
			elseif(ip1.eq.2) then
				aminoacid_name="NLYS"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NHIE"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			endif
		elseif(ip.eq.61) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="NASP"
			elseif(ip1.eq.3) then
				aminoacid_name="NCYT"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTYR"
			elseif(ip1.eq.7) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CARG"
			elseif(ip1.eq.2) then
				aminoacid_name="CLYS"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CHIE"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="CALA"
			endif
		elseif(ip.eq.62) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*3-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="CASP"
			elseif(ip1.eq.3) then
				aminoacid_name="CCYT"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 50

	elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
			group(ic)%gtype=="HIE") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="TYX") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
			group(ic)%gtype=="NHIE") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
			group(ic)%gtype=="CHIE") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
			ip=62
		endif

60		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="ARG"
			elseif(ip1.eq.2) then
				aminoacid_name="LYS"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="HIE"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			endif
		elseif(ip.eq.60) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="GLU"
			elseif(ip1.eq.2) then
				aminoacid_name="ASP"
			elseif(ip1.eq.3) then
				aminoacid_name="TYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CYT"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NARG"
			elseif(ip1.eq.2) then
				aminoacid_name="NLYS"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NHIE"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			endif
		elseif(ip.eq.61) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="NASP"
			elseif(ip1.eq.3) then
				aminoacid_name="NTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="NCYT"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CARG"
			elseif(ip1.eq.2) then
				aminoacid_name="CLYS"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*5-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CHIE"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip.eq.2) then
				aminoacid_name="CALA"
			endif
		elseif(ip.eq.62) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="CASP"
			elseif(ip1.eq.3) then
				aminoacid_name="CTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CCYT"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 60

	elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="ARG") then
			ip=30
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
			group(ic)%gtype=="HIE".or.group(ic)%gtype=="LYN") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="TYX") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NARG") then
			ip=31
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
			group(ic)%gtype=="NHIE".or.group(ic)%gtype=="NLYN") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CARG") then
			ip=32
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
			group(ic)%gtype=="CHIE".or.group(ic)%gtype=="CLYN") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
			ip=62
		endif

70		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.30) then
			aminoacid_name="ARG"
			goto 100
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="HIE"
			elseif(ip1.eq.6) then
				aminoacid_name="LYN"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			endif
		elseif(ip.eq.60) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="GLU"
			elseif(ip1.eq.2) then
				aminoacid_name="ASP"
			elseif(ip1.eq.3) then
				aminoacid_name="TYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CYT"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.31) then
			aminoacid_name="NARG"
			goto 100
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NHIE"
			elseif(ip1.eq.6) then
				aminoacid_name="NLYN"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			endif
		elseif(ip.eq.61) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="NASP"
			elseif(ip1.eq.3) then
				aminoacid_name="NTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="NCYT"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.32) then
			aminoacid_name="CARG"
			goto 100
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CHIE"
			elseif(ip1.eq.6) then
				aminoacid_name="CLYN"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip.eq.2) then
				aminoacid_name="CALA"
			endif
		elseif(ip.eq.62) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="CASP"
			elseif(ip1.eq.3) then
				aminoacid_name="CTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CCYT"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 70

	elseif(ph_value.ge.12.5) then
		if(group(ic)%gtype=="GLY") then
			ip=10
		elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or.group(ic)%gtype=="MET".or. &
			group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP") then
			ip=20
		elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or.group(ic)%gtype=="GLN".or. &
			group(ic)%gtype=="HIE".or.group(ic)%gtype=="LYN".or.group(ic)%gtype=="ARN") then
			ip=40
		elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA") then
			ip=50
		elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or.group(ic)%gtype=="TYX") then
			ip=60
		elseif(group(ic)%gtype=="NGLY") then
			ip=11
		elseif(group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or. &
			group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP") then
			ip=21
		elseif(group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. &
			group(ic)%gtype=="NHIE".or.group(ic)%gtype=="NLYN".or.group(ic)%gtype=="NARN") then
			ip=41
		elseif(group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NALA") then
			ip=51
		elseif(group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX") then
			ip=61
		elseif(group(ic)%gtype=="CGLY") then
			ip=12
		elseif(group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET".or. &
			group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
			ip=22
		elseif(group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. &
			group(ic)%gtype=="CHIE".or.group(ic)%gtype=="CLYN".or.group(ic)%gtype=="CARN") then
			ip=42
		elseif(group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
			ip=52
		elseif(group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
			ip=62
		endif

80		continue
		if(ip.eq.10) then
			aminoacid_name="GLY"
			goto 100
		elseif(ip.eq.20) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="LEU"
			elseif(ip1.eq.2) then
				aminoacid_name="VAL"
			elseif(ip1.eq.3) then
				aminoacid_name="ILE"
			elseif(ip1.eq.4) then
				aminoacid_name="MET"
			elseif(ip1.eq.5) then
				aminoacid_name="PHE"
			elseif(ip1.eq.6) then
				aminoacid_name="TRP"
			endif
		elseif(ip.eq.40) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="SER"
			elseif(ip1.eq.2) then
				aminoacid_name="THR"
			elseif(ip1.eq.3) then
				aminoacid_name="ASN"
			elseif(ip1.eq.4) then
				aminoacid_name="GLN"
			elseif(ip1.eq.5) then
				aminoacid_name="HIE"
			elseif(ip1.eq.6) then
				aminoacid_name="LYN"
			elseif(ip1.eq.7) then
				aminoacid_name="ARN"
			endif
		elseif(ip.eq.50) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="PRO"
			elseif(ip1.eq.2) then
				aminoacid_name="ALA"
			endif
		elseif(ip.eq.60) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="GLU"
			elseif(ip1.eq.2) then
				aminoacid_name="ASP"
			elseif(ip1.eq.3) then
				aminoacid_name="TYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CYT"
			endif
		elseif(ip.eq.11) then
			aminoacid_name="NGLY"
			goto 100
		elseif(ip.eq.21) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="NVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="NILE"
			elseif(ip1.eq.4) then
				aminoacid_name="NMET"
			elseif(ip1.eq.5) then
				aminoacid_name="NPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="NTRP"
			endif
		elseif(ip.eq.41) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NSER"
			elseif(ip1.eq.2) then
				aminoacid_name="NTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="NASN"
			elseif(ip1.eq.4) then
				aminoacid_name="NGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="NHIE"
			elseif(ip1.eq.6) then
				aminoacid_name="NLYN"
			elseif(ip1.eq.7) then
				aminoacid_name="NARN"
			endif
		elseif(ip.eq.51) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NPRO"
			elseif(ip1.eq.2) then
				aminoacid_name="NALA"
			endif
		elseif(ip.eq.61) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="NGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="NASP"
			elseif(ip1.eq.3) then
				aminoacid_name="NTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="NCYT"
			endif
		elseif(ip.eq.12) then
			aminoacid_name="CGLY"
			goto 100
		elseif(ip.eq.22) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*6-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CLEU"
			elseif(ip1.eq.2) then
				aminoacid_name="CVAL"
			elseif(ip1.eq.3) then
				aminoacid_name="CILE"
			elseif(ip1.eq.4) then
				aminoacid_name="CMET"
			elseif(ip1.eq.5) then
				aminoacid_name="CPHE"
			elseif(ip1.eq.6) then
				aminoacid_name="CTRP"
			endif
		elseif(ip.eq.42) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*7-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CSER"
			elseif(ip1.eq.2) then
				aminoacid_name="CTHR"
			elseif(ip1.eq.3) then
				aminoacid_name="CASN"
			elseif(ip1.eq.4) then
				aminoacid_name="CGLN"
			elseif(ip1.eq.5) then
				aminoacid_name="CHIE"
			elseif(ip1.eq.6) then
				aminoacid_name="CLYN"
			elseif(ip1.eq.7) then
				aminoacid_name="CARN"
			endif
		elseif(ip.eq.52) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*2-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CPRO"
			elseif(ip.eq.2) then
				aminoacid_name="CALA"
			endif
		elseif(ip.eq.62) then
			call ran_gen(ran2,zero)
			ip1=int(ran2*4-1.0e-3)+1
			if(ip1.eq.1) then
				aminoacid_name="CGLU"
			elseif(ip1.eq.2) then
				aminoacid_name="CASP"
			elseif(ip1.eq.3) then
				aminoacid_name="CTYX"
			elseif(ip1.eq.4) then
				aminoacid_name="CCYT"
			endif
		endif
		if(aminoacid_name==group(ic)%gtype) goto 80

	endif

100	continue

	return
	end subroutine mc_choose_aminoacid

	subroutine groupinfo(name, group_name, flag)
	implicit none
	integer*8					:: i, flag
	character*4				:: name, group_name(3)

	if(name=="GLY".or.name=="NGLY".or.name=="CGLY") then
		group_name(1)="GLY"
		group_name(2)="NGLY"
		group_name(3)="CGLY"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LEU".or.name=="NLEU".or.name=="CLEU") then
		group_name(1)="LEU"
		group_name(2)="NLEU"
		group_name(3)="CLEU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="VAL".or.name=="NVAL".or.name=="CVAL") then
		group_name(1)="VAL"
		group_name(2)="NVAL"
		group_name(3)="CVAL"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ILE".or.name=="NILE".or.name=="CILE") then
		group_name(1)="ILE"
		group_name(2)="NILE"
		group_name(3)="CILE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="MET".or.name=="NMET".or.name=="CMET") then
		group_name(1)="MET"
		group_name(2)="NMET"
		group_name(3)="CMET"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PHE".or.name=="NPHE".or.name=="CPHE") then
		group_name(1)="PHE"
		group_name(2)="NPHE"
		group_name(3)="CPHE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYR".or.name=="NTYR".or.name=="CTYR") then
		group_name(1)="TYR"
		group_name(2)="NTYR"
		group_name(3)="CTYR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYX".or.name=="NTYX".or.name=="CTYX") then
		group_name(1)="TYX"
		group_name(2)="NTYX"
		group_name(3)="CTYX"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TRP".or.name=="NTRP".or.name=="CTRP") then
		group_name(1)="TRP"
		group_name(2)="NTRP"
		group_name(3)="CTRP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARG".or.name=="NARG".or.name=="CARG") then
		group_name(1)="ARG"
		group_name(2)="NARG"
		group_name(3)="CARG"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARN".or.name=="NARN".or.name=="CARN") then
		group_name(1)="ARN"
		group_name(2)="NARN"
		group_name(3)="CARN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYN".or.name=="NLYN".or.name=="CLYN") then
		group_name(1)="LYN"
		group_name(2)="NLYN"
		group_name(3)="CLYN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYS".or.name=="NLYS".or.name=="CLYS") then
		group_name(1)="LYS"
		group_name(2)="NLYS"
		group_name(3)="CLYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="SER".or.name=="NSER".or.name=="CSER") then
		group_name(1)="SER"
		group_name(2)="NSER"
		group_name(3)="CSER"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="THR".or.name=="NTHR".or.name=="CTHR") then
		group_name(1)="THR"
		group_name(2)="NTHR"
		group_name(3)="CTHR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASN".or.name=="NASN".or.name=="CASN") then
		group_name(1)="ASN"
		group_name(2)="NASN"
		group_name(3)="CASN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLN".or.name=="NGLN".or.name=="CGLN") then
		group_name(1)="GLN"
		group_name(2)="NGLN"
		group_name(3)="CGLN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIE".or.name=="NHIE".or.name=="CHIE") then
		group_name(1)="HIE"
		group_name(2)="NHIE"
		group_name(3)="CHIE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIP".or.name=="NHIP".or.name=="CHIP") then
		group_name(1)="HIP"
		group_name(2)="NHIP"
		group_name(3)="CHIP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
		group_name(1)="PRO"
		group_name(2)="NPRO"
		group_name(3)="CPRO"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYS".or.name=="NCYS".or.name=="CCYS") then
		group_name(1)="CYS"
		group_name(2)="NCYS"
		group_name(3)="CCYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYT".or.name=="NCYT".or.name=="CCYT") then
		group_name(1)="CYT"
		group_name(2)="NCYT"
		group_name(3)="CCYT"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ALA".or.name=="NALA".or.name=="CALA") then
		group_name(1)="ALA"
		group_name(2)="NALA"
		group_name(3)="CALA"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLH".or.name=="NGLH".or.name=="CGLH") then
		group_name(1)="GLH"
		group_name(2)="NGLH"
		group_name(3)="CGLH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLU".or.name=="NGLU".or.name=="CGLU") then
		group_name(1)="GLU"
		group_name(2)="NGLU"
		group_name(3)="CGLU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASH".or.name=="NASH".or.name=="CASH") then
		group_name(1)="ASH"
		group_name(2)="NASH"
		group_name(3)="CASH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASP".or.name=="NASP".or.name=="CASP") then
		group_name(1)="ASP"
		group_name(2)="NASP"
		group_name(3)="CASP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	endif
5	continue

	return
	end subroutine groupinfo

	subroutine scmf_choose_aminoacid(ip, aminoacid_number, aminoacid_name)
	implicit none
	integer*8					:: aminoacid_number
	integer*8					:: ip, ip1, ip2, i
	double precision					:: ran2
	character*4				:: char, aminoacid_name(10)

	if(ph_value.le.3.9) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=7
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TYR"
			aminoacid_name(7)="TRP"
		elseif(ip==3) then
			aminoacid_number=3
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
			aminoacid_name(3)="HIP"
		elseif(ip==4) then
			aminoacid_number=6
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="GLH"
			aminoacid_name(6)="ASH"
		elseif(ip==5) then
			aminoacid_number=3
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
			aminoacid_name(3)="CYS"
		endif

	elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=7
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TYR"
			aminoacid_name(7)="TRP"
		elseif(ip==3) then
			aminoacid_number=3
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
			aminoacid_name(3)="HIP"
		elseif(ip==4) then
			aminoacid_number=5
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="GLH"
		elseif(ip==5) then
			aminoacid_number=3
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
			aminoacid_name(3)="CYS"
		elseif(ip==6) then
			aminoacid_number=1
			aminoacid_name(1)="ASP"
		endif

	elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=7
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TYR"
			aminoacid_name(7)="TRP"
		elseif(ip==3) then
			aminoacid_number=3
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
			aminoacid_name(3)="HIP"
		elseif(ip==4) then
			aminoacid_number=4
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
		elseif(ip==5) then
!			aminoacid_number=3
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
!			aminoacid_name(3)="CYS"
		elseif(ip==6) then
			aminoacid_number=2
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
		endif

	elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=7
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TYR"
			aminoacid_name(7)="TRP"
		elseif(ip==3) then
			aminoacid_number=2
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
		elseif(ip==4) then
			aminoacid_number=5
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="HIE"
		elseif(ip==5) then
!			aminoacid_number=3
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
!			aminoacid_name(3)="CYS"
		elseif(ip==6) then
			aminoacid_number=2
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
		endif

	elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=7
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TYR"
			aminoacid_name(7)="TRP"
		elseif(ip==3) then
			aminoacid_number=2
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
		elseif(ip==4) then
			aminoacid_number=5
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="HIE"
		elseif(ip==5) then
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
		elseif(ip==6) then
			aminoacid_number=3
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
			aminoacid_name(3)="CYT"
		endif

	elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=6
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TRP"
		elseif(ip==3) then
			aminoacid_number=2
			aminoacid_name(1)="ARG"
			aminoacid_name(2)="LYS"
		elseif(ip==4) then
			aminoacid_number=5
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="HIE"
		elseif(ip==5) then
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
		elseif(ip==6) then
			aminoacid_number=4
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
			aminoacid_name(3)="TYX"
			aminoacid_name(4)="CYT"
		endif

	elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=6
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TRP"
		elseif(ip==3) then
			aminoacid_number=1
			aminoacid_name(1)="ARG"
		elseif(ip==4) then
			aminoacid_number=6
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="HIE"
			aminoacid_name(6)="LYN"
		elseif(ip==5) then
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
		elseif(ip==6) then
			aminoacid_number=4
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
			aminoacid_name(3)="TYX"
			aminoacid_name(4)="CYT"
		endif

	elseif(ph_value.ge.12.5) then
		if(ip==1) then
			aminoacid_number=1
			aminoacid_name(1)="GLY"
		elseif(ip==2) then
			aminoacid_number=6
			aminoacid_name(1)="LEU"
			aminoacid_name(2)="VAL"
			aminoacid_name(3)="ILE"
			aminoacid_name(4)="MET"
			aminoacid_name(5)="PHE"
			aminoacid_name(6)="TRP"
		elseif(ip==4) then
			aminoacid_number=7
			aminoacid_name(1)="SER"
			aminoacid_name(2)="THR"
			aminoacid_name(3)="ASN"
			aminoacid_name(4)="GLN"
			aminoacid_name(5)="HIE"
			aminoacid_name(6)="LYN"
			aminoacid_name(7)="ARN"
		elseif(ip==5) then
			aminoacid_number=2
			aminoacid_name(1)="PRO"
			aminoacid_name(2)="ALA"
		elseif(ip==6) then
			aminoacid_number=4
			aminoacid_name(1)="GLU"
			aminoacid_name(2)="ASP"
			aminoacid_name(3)="TYX"
			aminoacid_name(4)="CYT"
		endif

	endif

	do i=1, (aminoacid_number-1)
		call ran_gen(ran2,zero)
		ip1=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip1.gt.aminoacid_number) ip1=aminoacid_number
		call ran_gen(ran2,zero)
		ip2=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip2.gt.aminoacid_number) ip2=aminoacid_number

		char=aminoacid_name(ip2)
		aminoacid_name(ip2)=aminoacid_name(ip1)
		aminoacid_name(ip1)=char
	enddo

	return
	end subroutine scmf_choose_aminoacid

	subroutine sidechain_category(ic, group, Iclass, grade, grade_num, index, monitor)
	implicit none
	integer*8								:: grade, grade_num(6), monitor(6), i, ic
	type(groupdetails)					:: group(gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6)

	grade_num=0
	if(group(ic)%gtype=="VAL".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="CVAL") then
		grade=1
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1)
				Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2)
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo
	elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="CLEU") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="ILE".or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="CILE") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB".or.group(ic)%atype2(i)=="CG2".or.group(ic)%atype2(i)=="HG21" &
				.or.group(ic)%atype2(i)=="HG22".or.group(ic)%atype2(i)=="HG23".or.group(ic)%atype2(i)=="CG1") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG1") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="PHE".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="CPHE") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="TRP".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CTRP") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="TYR".or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="CTYR") then
		grade=3
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HH") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="OH") monitor(3)=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="TYX".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CTYX") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="OH") monitor(3)=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="CSER") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="OG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="OG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="THR".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="CTHR") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HG1") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="OG1") monitor(2)=grade_num(2)
			endif
		enddo
	elseif(group(ic)%gtype=="CYS".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="CCYS") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="SG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="SG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="CYT".or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="CCYT") then
		grade=1
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="SG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="SG") monitor(2)=grade_num(2)
			endif
		enddo
	elseif(group(ic)%gtype=="MET".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="CMET") then
		grade=3
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="SD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="SD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ic)%gtype=="ASN".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="CASN") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="GLN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CGLN") then
		grade=3
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="ASH".or.group(ic)%gtype=="NASH".or.group(ic)%gtype=="CASH") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="CGLU") then
		grade=3
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ic)%gtype=="GLH".or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="CGLH") then
		grade=3
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ic)%gtype=="HIE".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CHIE") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="HIP".or.group(ic)%gtype=="NHIP".or.group(ic)%gtype=="CHIP") then
		grade=2
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ic)%gtype=="LYS".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CLYS") then
		grade=4
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ic)%atype2(i)=="CE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ic)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=group(ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(group(ic)%gtype=="LYN".or.group(ic)%gtype=="NLYN".or.group(ic)%gtype=="CLYN") then
		grade=4
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ic)%atype2(i)=="CE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ic)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=group(ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG") then
		grade=4
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="NE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ic)%atype2(i)=="NE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ic)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=group(ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(group(ic)%gtype=="ARN".or.group(ic)%gtype=="NARN".or.group(ic)%gtype=="CARN") then
		grade=4
		do i=1, group(ic)%cnum2
			if(group(ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ic)%coo2(i,2); 
				Iclass(1)%member(grade_num(1),3)=group(ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ic)%atype2(i)=="HB2".or.group(ic)%atype2(i)=="HB3".or.group(ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ic)%coo2(i,2); 
				Iclass(2)%member(grade_num(2),3)=group(ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ic)%atype2(i)=="HG2".or.group(ic)%atype2(i)=="HG3".or.group(ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ic)%coo2(i,2); 
				Iclass(3)%member(grade_num(3),3)=group(ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ic)%atype2(i)=="HD2".or.group(ic)%atype2(i)=="HD3".or.group(ic)%atype2(i)=="NE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ic)%coo2(i,2); 
				Iclass(4)%member(grade_num(4),3)=group(ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ic)%atype2(i)=="NE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ic)%coo2(i,2); 
				Iclass(5)%member(grade_num(5),3)=group(ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	endif

	return
	end subroutine sidechain_category

	subroutine dihedralangle_reading(gtype, dihedral_num, dihedral)
	implicit none
	integer*8								:: dihedral_num, i, j
	character*4							:: gtype
	type(dihedralparameters)			:: dihedral

	open(10, file='~/PepBD/lib/DihedralAngle/'//trim(gtype), status="old")
		read(10, "(i8)") dihedral_num
		do i=1, dihedral_num
			read(10,"(5i8)") dihedral%iph(i), dihedral%jph(i), dihedral%kph(i), dihedral%lph(i), dihedral%multiply(i)
			do j=1, dihedral%multiply(i)
				read(10,"(3e16.8)") dihedral%pk(i,j), dihedral%pn(i,j), dihedral%phase(i,j)
			enddo
		enddo
	close(10)

	return
	end subroutine dihedralangle_reading

end module database

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module transplant

	use constant
	use mathfunction

	contains
	subroutine residue_replace(ic, group, groupdata_backup, ip, aa_group, temp_group)
	implicit none
	integer*8								:: ic, ip, i, j, flag
	type(groupdetails)					:: group(gnum), temp_group(gnum), aa_group(40)
	type(databackup)					:: groupdata_backup(gnum)

	temp_group=group
	if(temp_group(ic)%gtype=="PRO".or.temp_group(ic)%gtype=="NPRO".or.temp_group(ic)%gtype=="CPRO") then
		if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			flag=0
			do i=1, temp_group(ic)%cnum1
				if(temp_group(ic)%atype1(i)=="H2".or.temp_group(ic)%atype1(i)=="H3") then
					flag=1
				endif
			enddo

			do i=1, aa_group(ip)%cnum1
				if(aa_group(ip)%atype1(i)=="H") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="N") then
							if(flag==1) then
								temp_group(ic)%atype1(j+1)="H1"
							else
								temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							endif
							temp_group(ic)%coo1((j+1),1)=groupdata_backup(ic)%coo(1)
							temp_group(ic)%coo1((j+1),2)=groupdata_backup(ic)%coo(2)
							temp_group(ic)%coo1((j+1),3)=groupdata_backup(ic)%coo(3)
							goto 10
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				elseif(aa_group(ip)%atype1(i)=="HA2") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="CA") then
							temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							temp_group(ic)%coo1((j+1),1)=aa_group(ip)%coo1(i,1)
							temp_group(ic)%coo1((j+1),2)=aa_group(ip)%coo1(i,2)
							temp_group(ic)%coo1((j+1),3)=aa_group(ip)%coo1(i,3)
							goto 10
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				elseif(aa_group(ip)%atype1(i)=="HA3") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="HA2") then
							temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							temp_group(ic)%coo1((j+1),1)=aa_group(ip)%coo1(i,1)
							temp_group(ic)%coo1((j+1),2)=aa_group(ip)%coo1(i,2)
							temp_group(ic)%coo1((j+1),3)=aa_group(ip)%coo1(i,3)
							goto 10
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				endif
10				continue
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
		else
			temp_group(ic)%gtype=aa_group(ip)%gtype
			flag=0
			do i=1, temp_group(ic)%cnum1
				if(temp_group(ic)%atype1(i)=="H2".or.temp_group(ic)%atype1(i)=="H3") then
					flag=1
				endif
			enddo

			do i=1, aa_group(ip)%cnum1
				if(aa_group(ip)%atype1(i)=="H") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="N") then
							if(flag==1) then
								temp_group(ic)%atype1(j+1)="H1"
							else
								temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							endif
							temp_group(ic)%coo1((j+1),1)=groupdata_backup(ic)%coo(1)
							temp_group(ic)%coo1((j+1),2)=groupdata_backup(ic)%coo(2)
							temp_group(ic)%coo1((j+1),3)=groupdata_backup(ic)%coo(3)
							goto 20
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				endif
			enddo
20			continue
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		endif
	elseif(temp_group(ic)%gtype=="GLY".or.temp_group(ic)%gtype=="NGLY".or.temp_group(ic)%gtype=="CGLY") then
		if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			flag=0
			do i=1, (temp_group(ic)%cnum1-1)
				if(temp_group(ic)%atype1(i)=="H1".or.temp_group(ic)%atype1(i)=="H".or.flag==1) then
					temp_group(ic)%atype1(i)=temp_group(ic)%atype1(i+1)
					temp_group(ic)%coo1(i,1)=temp_group(ic)%coo1((i+1),1)
					temp_group(ic)%coo1(i,2)=temp_group(ic)%coo1((i+1),2)
					temp_group(ic)%coo1(i,3)=temp_group(ic)%coo1((i+1),3)
					flag=1
				endif
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1

			flag=0
			do i=1, (temp_group(ic)%cnum1-1)
				if(temp_group(ic)%atype1(i)=="HA2") then
					temp_group(ic)%atype1(i)="HA"
				elseif(temp_group(ic)%atype1(i)=="HA3".or.flag==1) then
					temp_group(ic)%atype1(i)=temp_group(ic)%atype1(i+1)
					temp_group(ic)%coo1(i,1)=temp_group(ic)%coo1((i+1),1)
					temp_group(ic)%coo1(i,2)=temp_group(ic)%coo1((i+1),2)
					temp_group(ic)%coo1(i,3)=temp_group(ic)%coo1((i+1),3)
					flag=1
				endif
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		else
			temp_group(ic)%gtype=aa_group(ip)%gtype
			flag=0
			do i=1, (temp_group(ic)%cnum1-1)
				if(temp_group(ic)%atype1(i)=="HA2") then
					temp_group(ic)%atype1(i)="HA"
				elseif(temp_group(ic)%atype1(i)=="HA3".or.flag==1) then
					temp_group(ic)%atype1(i)=temp_group(ic)%atype1(i+1)
					temp_group(ic)%coo1(i,1)=temp_group(ic)%coo1((i+1),1)
					temp_group(ic)%coo1(i,2)=temp_group(ic)%coo1((i+1),2)
					temp_group(ic)%coo1(i,3)=temp_group(ic)%coo1((i+1),3)
					flag=1
				endif
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		endif
	else
		if(aa_group(ip)%gtype=="PRO".or.aa_group(ip)%gtype=="NPRO".or.aa_group(ip)%gtype=="CPRO") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			flag=0
			do i=1, (temp_group(ic)%cnum1-1)
				if(temp_group(ic)%atype1(i)=="H1".or.temp_group(ic)%atype1(i)=="H".or.flag==1) then
					temp_group(ic)%atype1(i)=temp_group(ic)%atype1(i+1)
					temp_group(ic)%coo1(i,1)=temp_group(ic)%coo1((i+1),1)
					temp_group(ic)%coo1(i,2)=temp_group(ic)%coo1((i+1),2)
					temp_group(ic)%coo1(i,3)=temp_group(ic)%coo1((i+1),3)
					flag=1
				endif
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		elseif(aa_group(ip)%gtype=="GLY".or.aa_group(ip)%gtype=="NGLY".or.aa_group(ip)%gtype=="CGLY") then
			temp_group(ic)%gtype=aa_group(ip)%gtype
			do i=1, aa_group(ip)%cnum1
				if(aa_group(ip)%atype1(i)=="HA2") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="CA") then
							temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							temp_group(ic)%coo1((j+1),1)=aa_group(ip)%coo1(i,1)
							temp_group(ic)%coo1((j+1),2)=aa_group(ip)%coo1(i,2)
							temp_group(ic)%coo1((j+1),3)=aa_group(ip)%coo1(i,3)
							goto 30
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				elseif(aa_group(ip)%atype1(i)=="HA3") then
					temp_group(ic)%cnum1=temp_group(ic)%cnum1+1
					do j=(temp_group(ic)%cnum1-1), 1, -1
						if(temp_group(ic)%atype1(j)=="HA2") then
							temp_group(ic)%atype1(j+1)=aa_group(ip)%atype1(i)
							temp_group(ic)%coo1((j+1),1)=aa_group(ip)%coo1(i,1)
							temp_group(ic)%coo1((j+1),2)=aa_group(ip)%coo1(i,2)
							temp_group(ic)%coo1((j+1),3)=aa_group(ip)%coo1(i,3)
							goto 30
						else
							temp_group(ic)%atype1(j+1)=temp_group(ic)%atype1(j)
							temp_group(ic)%coo1((j+1),1)=temp_group(ic)%coo1(j,1)
							temp_group(ic)%coo1((j+1),2)=temp_group(ic)%coo1(j,2)
							temp_group(ic)%coo1((j+1),3)=temp_group(ic)%coo1(j,3)
						endif
					enddo
				endif
30				continue
			enddo
			temp_group(ic)%cnum1=temp_group(ic)%cnum1-1
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
		else
			temp_group(ic)%gtype=aa_group(ip)%gtype
			temp_group(ic)%cnum2=aa_group(ip)%cnum2
			do i=1, aa_group(ip)%cnum2
				temp_group(ic)%atype2(i)=aa_group(ip)%atype2(i)
				temp_group(ic)%coo2(i,1)=aa_group(ip)%coo2(i,1)
				temp_group(ic)%coo2(i,2)=aa_group(ip)%coo2(i,2)
				temp_group(ic)%coo2(i,3)=aa_group(ip)%coo2(i,3)
			enddo
		endif
	endif

	return
	end subroutine residue_replace

	subroutine check_transplant(flag, num1, num2, group, feedback)
	implicit none
	integer*8							:: i, j, k, num1, num2, ic1, ic2, flag, feedback
	double precision							:: rij
	type(groupdetails)				:: group(gnum)

	feedback=1
	if(flag==0) then
		ic1=num1
		do i=1, group(ic1)%cnum2
			do k=1, gnum
				if (k.eq.ic1) goto 10
				do j=1, group(k)%cnum1
					if(k==(ic1+1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype1(j)=="N") goto 30
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
							goto 30
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
					rij=sqrt(rij)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(rij<1.55) then
							feedback=0
							goto 20
						endif
					elseif(rij<1.55) then
						feedback=0
						goto 20
					endif
30					continue
				enddo
				do j=1, group(k)%cnum2
					rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
					rij=sqrt(rij)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(rij<1.55) then
							feedback=0
							goto 20
						endif
					elseif(rij<1.55) then
						feedback=0
						goto 20
					endif
				enddo
				do j=1, group(k)%cnum3
					if(k==(ic1-1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype3(j)=="C") goto 40
					if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
						if(group(ic1)%atype2(i)=="CD") then
							goto 40
						elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
							if(group(k)%atype3(j)=="C") goto 40
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
					rij=sqrt(rij)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(rij<1.55) then
							feedback=0
							goto 20
						endif
					elseif(rij<1.55) then
						feedback=0
						goto 20
					endif
40				    continue
				enddo
10				continue
			enddo
		enddo
20		continue
	elseif(flag==1) then
		ic1=num1
		ic2=num2
		do i=1, group(ic1)%cnum2
			do k=1, gnum
				if (k.eq.ic1) goto 50
				if (k.eq.ic2) then
					do j=1, group(k)%cnum1
						if(k==(ic1+1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype1(j)=="N") goto 53
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
								goto 53
							endif
						endif
						rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
						    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 60
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 60
						endif
53						continue
					enddo
					do j=1, group(k)%cnum3
						if(k==(ic1-1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype3(j)=="C") goto 57
						if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
							if(group(ic1)%atype2(i)=="CD") then
								goto 57
							elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
								if(group(k)%atype3(j)=="C") goto 57
							endif
						endif
						rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
						    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 60
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 60
						endif
57					    continue
					enddo
				else
					do j=1, group(k)%cnum1
						if(k==(ic1+1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype1(j)=="N") goto 70
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
								goto 70
							endif
						endif
						rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
						    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 60
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 60
						endif
70						continue
					enddo
					do j=1, group(k)%cnum2
						rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
						    (group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 60
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 60
						endif
					enddo
					do j=1, group(k)%cnum3
						if(k==(ic1-1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype3(j)=="C") goto 80
						if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
							if(group(ic1)%atype2(i)=="CD") then
								goto 80
							elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
								if(group(k)%atype3(j)=="C") goto 80
							endif
						endif
						rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
						    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 60
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 60
						endif
80					    continue
					enddo
				endif
50				continue
			enddo
		enddo
60		continue
	elseif(flag==2) then
		do ic1=1, sitenum
			do i=1, group(ic1)%cnum1
				do k=ic1+1, gnum
					do j=1, group(k)%cnum2
						rij=(group(ic1)%coo1(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic1)%coo1(i,3)-group(k)%coo2(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
					enddo

					do j=1, helix_num
						if(ic1.ge.helix_start(j).and.ic1.le.helix_end(j).and.k.ge.helix_start(j).and.k.le.helix_end(j)) goto 97
					enddo

					do j=1, group(k)%cnum1
						if(k==(ic1+1)) goto 95
						rij=(group(ic1)%coo1(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic1)%coo1(i,3)-group(k)%coo1(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
95						continue
					enddo
					do j=1, group(k)%cnum3
						rij=(group(ic1)%coo1(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic1)%coo1(i,3)-group(k)%coo3(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
					enddo
97					continue
				enddo
			enddo

			do i=1, group(ic1)%cnum2
				do k=ic1+1, gnum
					do j=1, group(k)%cnum1
						if(k==(ic1+1).and.group(ic1)%atype2(i)=="CB".and.group(k)%atype1(j)=="N") goto 100
						rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
100						continue
					enddo
					do j=1, group(k)%cnum2
						rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
						   group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
					enddo
					do j=1, group(k)%cnum3
						rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
						rij=sqrt(rij)
						if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
					enddo
				enddo
			enddo

			do i=1, group(ic1)%cnum3
				do k=ic1+1, gnum
					do j=1, group(k)%cnum2
						if(k==(ic1+1)) goto 110
						rij=(group(ic1)%coo3(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic1)%coo3(i,3)-group(k)%coo2(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
110						continue
					enddo

					do j=1, helix_num
						if(ic1.ge.helix_start(j).and.ic1.le.helix_end(j).and.k.ge.helix_start(j).and.k.le.helix_end(j)) goto 120
					enddo

					do j=1, group(k)%cnum1
						if(k==(ic1+1)) goto 105
						rij=(group(ic1)%coo3(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic1)%coo3(i,3)-group(k)%coo1(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
105						continue
					enddo
					do j=1, group(k)%cnum3
						if(k==(ic1+1)) goto 115
						rij=(group(ic1)%coo3(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic1)%coo3(i,3)-group(k)%coo3(j,3))**2
						rij=sqrt(rij)
						if(group(k)%gtype=="PRO".or.group(k)%gtype=="NPRO".or.group(k)%gtype=="CPRO") then
							if(rij<1.55) then
								feedback=0
								goto 90
							endif
						elseif(rij<1.55) then
							feedback=0
							goto 90
						endif
115						continue
					enddo
120					continue
				enddo
			enddo

		enddo
90		continue
	endif

	return
	end subroutine 	check_transplant

	subroutine check_backbone(group, feedback)
	!this ensures there are no steric clashes between peptide backbone with either itself or the receptor
	implicit none
	integer*8				:: feedback
	integer*8				:: i, j, k, ic1
	double precision					:: rij
	real					:: min_rij = 1.6*1.6
	type(groupdetails)		:: group(gnum)

	feedback=1
	do ic1=1, sitenum
		do i=1, group(ic1)%cnum1
			do k=ic1+1, gnum
				do j=1, helix_num
					if(ic1.ge.helix_start(j).and.ic1.le.helix_end(j).and.k.ge.helix_start(j).and.k.le.helix_end(j)) goto 20
				enddo

				if(k.gt.(ic1+1).or.k.gt.sitenum) then
					do j=1, group(k)%cnum1
						rij=(group(ic1)%coo1(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic1)%coo1(i,3)-group(k)%coo1(j,3))**2
						if(rij<min_rij) then
							feedback=0
							goto 10
						endif
					enddo

					do j=1, group(k)%cnum3
						rij=(group(ic1)%coo1(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic1)%coo1(i,3)-group(k)%coo3(j,3))**2
						if(rij<min_rij) then
							feedback=0
							goto 10
						endif
					enddo

					if(k.gt.sitenum) then
						do j=1, group(k)%cnum2
							rij=(group(ic1)%coo1(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo1(i,2)-group(k)%coo2(j,2))**2+ &
								(group(ic1)%coo1(i,3)-group(k)%coo2(j,3))**2
							if(rij<min_rij) then
								feedback=0
								goto 10
							endif
						enddo
					endif
				endif
20				continue
			enddo
		enddo

		do i=1, group(ic1)%cnum3
			do k=ic1+1, gnum
				do j=1, helix_num
					if(ic1.ge.helix_start(j).and.ic1.le.helix_end(j).and.k.ge.helix_start(j).and.k.le.helix_end(j)) goto 30
				enddo

				if(k.gt.(ic1+1).or.k.gt.sitenum) then
					do j=1, group(k)%cnum1
						rij=(group(ic1)%coo3(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic1)%coo3(i,3)-group(k)%coo1(j,3))**2
						if(rij<min_rij) then
							feedback=0
							goto 10
						endif
					enddo

					do j=1, group(k)%cnum3
						rij=(group(ic1)%coo3(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic1)%coo3(i,3)-group(k)%coo3(j,3))**2
						if(rij<min_rij) then
							feedback=0
							goto 10
						endif
					enddo

					if(k.gt.sitenum) then
						do j=1, group(k)%cnum2
							rij=(group(ic1)%coo3(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo3(i,2)-group(k)%coo2(j,2))**2+ &
								(group(ic1)%coo3(i,3)-group(k)%coo2(j,3))**2
							if(rij<min_rij) then
								feedback=0
								goto 10
							endif
						enddo
					endif
				endif
30				continue
			enddo
		enddo
	enddo
10	continue

	return
	end subroutine check_backbone

end module transplant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module energy_calculation

	use constant
	use mathfunction
	use database

	contains
	subroutine vdwenergy_sidechain_receptor(flag, ic1, ic2, group, group_para, energy)
	implicit none
	integer*8							:: i, j, k, ic1, ic2, flag
	double precision							:: energy, rij, epsion, r0, acoeff, bcoeff, vdw
	type(groupdetails)				:: group(gnum)
	type(energyparameters)			:: group_para(gnum)

	energy=0.0

	if(flag==0) then
		do i=1, group(ic1)%cnum2
			do k=sitenum+1, gnum
				do j=1, group(k)%cnum1
					if(k==(ic1+1)) goto 20
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
							goto 20
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion1(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r1(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
20					continue
				enddo
				do j=1, group(k)%cnum2
					rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion2(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r2(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
				enddo
				do j=1, group(k)%cnum3
					if(k==(ic1-1)) goto 30
					if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
						if(group(ic1)%atype2(i)=="CD") then
							goto 30
						elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
							if(group(k)%atype3(j)=="C") goto 30
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion3(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r3(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
30					continue
				enddo
			enddo
		enddo
	elseif(flag==1) then
		do i=1, group(ic1)%cnum2
			do k=1, gnum
				if (k.eq.ic1.or.k.eq.ic2) goto 40
					do j=1, group(k)%cnum1
					if(k==(ic1+1)) goto 50
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
						if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
							goto 50
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion1(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r1(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
50					continue
				enddo
				do j=1, group(k)%cnum2
					rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion2(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r2(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
				enddo
				do j=1, group(k)%cnum3
					if(k==(ic1-1)) goto 60
					if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
						if(group(ic1)%atype2(i)=="CD") then
							goto 60
						elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
							if(group(k)%atype3(j)=="C") goto 60
						endif
					endif
					rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
					    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
					epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion3(j))
					r0=group_para(ic1)%r2(i)+group_para(k)%r3(j)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
					group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
						if(rij<10.0) rij=12.25
					endif
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
60					continue
				enddo
40				continue
			enddo
		enddo
	elseif(flag==2) then
		do i=1, group(ic1)%cnum2
			k=ic2
			do j=1, group(k)%cnum1
				if(k==(ic1+1)) goto 70
				if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO") then
					if(k==(ic1-1).and.group(ic1)%atype2(i)=="CD".and.group(k)%atype1(j)=="CA") then
						goto 70
					endif
				endif
				rij=(group(ic1)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
				    (group(ic1)%coo2(i,3)-group(k)%coo1(j,3))**2
				epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion1(j))
				r0=group_para(ic1)%r2(i)+group_para(k)%r1(j)
				acoeff=epsion*(r0**12)
				bcoeff=epsion*2*(r0**6)
				if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
				group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
					if(rij<10.0) rij=12.25
				endif
				vdw=acoeff/(rij**6)-bcoeff/(rij**3)
				energy=energy+vdw
70				continue
			enddo
			do j=1, group(k)%cnum2
				rij=(group(ic1)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
				    (group(ic1)%coo2(i,3)-group(k)%coo2(j,3))**2
				epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion2(j))
				r0=group_para(ic1)%r2(i)+group_para(k)%r2(j)
				acoeff=epsion*(r0**12)
				bcoeff=epsion*2*(r0**6)
				if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
				group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
					if(rij<10.0) rij=12.25
				endif
				vdw=acoeff/(rij**6)-bcoeff/(rij**3)
				energy=energy+vdw
			enddo
			do j=1, group(k)%cnum3
				if(k==(ic1-1)) goto 80
				if((group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO").and.k==(ic1-1)) then
					if(group(ic1)%atype2(i)=="CD") then
						goto 80
					elseif(group(ic1)%atype2(i)=="HD2".or.group(ic1)%atype2(i)=="HD3") then
						if(group(k)%atype3(j)=="C") goto 80
					endif
				endif
				rij=(group(ic1)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic1)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
				    (group(ic1)%coo2(i,3)-group(k)%coo3(j,3))**2
				epsion=sqrt(group_para(ic1)%epsion2(i)*group_para(k)%epsion3(j))
				r0=group_para(ic1)%r2(i)+group_para(k)%r3(j)
				acoeff=epsion*(r0**12)
				bcoeff=epsion*2*(r0**6)
				if(group(ic1)%gtype=="PRO".or.group(ic1)%gtype=="NPRO".or.group(ic1)%gtype=="CPRO".or. &
				group(ic1+1)%gtype=="PRO".or.group(ic1+1)%gtype=="NPRO".or.group(ic1+1)%gtype=="CPRO") then
					if(rij<10.0) rij=12.25
				endif
				vdw=acoeff/(rij**6)-bcoeff/(rij**3)
				energy=energy+vdw
80				continue
			enddo
		enddo
	endif

	return
	end subroutine vdwenergy_sidechain_receptor

	subroutine sidechain_energy(stage, ic, group, group_para, numex, inb, numex4, inb4, dihedral_num, dihedral, energy)
	implicit none
	integer*8								:: i, j, k, l, ic, i_id, j_id, flag, stage
	integer*8								:: numex(60), inb(60,20), numex4(60), inb4(60,60)
	integer*8								:: natom, dihedral_num, ip, jp, kp, lp
	double precision								:: energy, rij, epsion, r0, acoeff, bcoeff, dielecons4solute, vdw, ele
	double precision								:: p1(3), p2(3), p3(3), p4(3), angle, dihedral_energy, potential
	double precision							    :: mdcrd(atom_num,3)
	type(groupdetails)					:: group(gnum)
	type(energyparameters)				:: group_para(gnum)
	type(dihedralparameters)		    :: dihedral

	
	natom=0
	do i=1, group(ic)%cnum1
		natom=natom+1
		mdcrd(natom,1)=group(ic)%coo1(i,1)
		mdcrd(natom,2)=group(ic)%coo1(i,2)
		mdcrd(natom,3)=group(ic)%coo1(i,3)
	enddo
	do i=1, group(ic)%cnum2
		natom=natom+1
		mdcrd(natom,1)=group(ic)%coo2(i,1)
		mdcrd(natom,2)=group(ic)%coo2(i,2)
		mdcrd(natom,3)=group(ic)%coo2(i,3)
	enddo
	do i=1, group(ic)%cnum3
		natom=natom+1
		mdcrd(natom,1)=group(ic)%coo3(i,1)
		mdcrd(natom,2)=group(ic)%coo3(i,2)
		mdcrd(natom,3)=group(ic)%coo3(i,3)
	enddo

    !calculate dihedral energy for all dihedral angles in ligand
	!$OMP PARALLEL PRIVATE(p1, p2, p3, p4, angle, ip, jp, kp, lp, i, j, potential)
	dihedral_energy=0.0
    !$OMP DO REDUCTION (+ : dihedral_energy) 
	do i=1, dihedral_num
		ip=dihedral%iph(i); jp=dihedral%jph(i); kp=dihedral%kph(i); lp=dihedral%lph(i)
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		do j=1, dihedral%multiply(i)
			potential=dihedral%pk(i,j)*(1+cosd(dihedral%pn(i,j)*angle-dihedral%phase(i,j)))
			dihedral_energy=dihedral_energy+potential
		enddo
	enddo
    !$OMP END DO
	!$OMP END PARALLEL

    !get interaction energy between sidechain and rest of system
    !$OMP PARALLEL PRIVATE(i,j,k,l,i_id,j_id,flag,rij, epsion,r0,acoeff,bcoeff,vdw,ele,dielecons4solute)
	energy=0.0
	!$OMP DO REDUCTION (+ : energy) 
	do i=1, group(ic)%cnum2
		if(group(ic)%atype2(i)=="CB") goto 10
		do k=1, gnum
			if(stage==0) then
                !if stage == 0, also calculate energy of sidechain interact with its own residue, and only calculate vdW energy
				if(k.eq.ic) then
					i_id=group_para(ic)%atomid2(i)
					do j=1, group(k)%cnum1
						j_id=group_para(k)%atomid1(j)
						do l=1, numex(j_id)
							if(i_id.eq.inb(j_id,l)) goto 20
						enddo
						flag=0
						do l=1, numex4(j_id)
							if(i_id.eq.inb4(j_id,l)) then
								flag=1
								goto 30
							endif
						enddo
30						continue
						rij=(group(ic)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo1(j,3))**2
                        !probably could remove cutoff check here since these atoms are in same residue and should be within cutoff of each (perhaps except maybe for very long sidechains like Lysine? Could check at some point)
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion1(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r1(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(flag==1) then
                                vdw=vdw/vdw14_coeff
                            endif
                        endif
						energy = energy+vdw
20						continue
					enddo
					do j=1, group(k)%cnum2
						j_id=group_para(k)%atomid2(j)
						if(i_id.eq.j_id) goto 40
						do l=1, numex(j_id)
							if(i_id.eq.inb(j_id,l)) goto 40
						enddo
						flag=0
						do l=1, numex4(j_id)
							if(i_id.eq.inb4(j_id,l)) then
								flag=1
								goto 50
							endif
						enddo
50						continue
						rij=(group(ic)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo2(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion2(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r2(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(flag==1) then
                                vdw=vdw/vdw14_coeff
                            endif
                            energy =energy +vdw
                        endif
40						continue
					enddo
					do j=1, group(k)%cnum3
						j_id=group_para(k)%atomid3(j)
						do l=1, numex(j_id)
							if(i_id.eq.inb(j_id,l)) goto 60
						enddo
						flag=0
						do l=1, numex4(j_id)
							if(i_id.eq.inb4(j_id,l)) then
								flag=1
								goto 70
							endif
						enddo
70						continue
						rij=(group(ic)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo3(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion3(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r3(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(flag==1) then
                                vdw=vdw/vdw14_coeff
                            endif
                            energy =energy+vdw
                        endif
						
60						continue
					enddo
				else
					do j=1, group(k)%cnum1
						rij=(group(ic)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo1(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then 
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion1(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r1(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            energy =energy+vdw
                        endif
					enddo
					do j=1, group(k)%cnum2
						rij=(group(ic)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo2(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion2(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r2(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            energy =energy +vdw
                        endif
					enddo
					do j=1, group(k)%cnum3
						rij=(group(ic)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo3(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion3(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r3(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            energy =energy +vdw
                        endif
					enddo
				endif
			elseif(stage==1) then
                !only calculate interaction of sidechain with other residues (includes vdW and electrostatic energy)
				if(k.ne.ic) then
					do j=1, group(k)%cnum1
						rij=(group(ic)%coo2(i,1)-group(k)%coo1(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo1(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo1(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion1(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r1(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(group_para(ic)%dielecons2(i).ge.group_para(k)%dielecons1(j)) then
                                dielecons4solute=group_para(ic)%dielecons2(i)
                            else
                                dielecons4solute=group_para(k)%dielecons1(j)
                            endif
                            ele=(group_para(ic)%charge2(i)*group_para(k)%charge1(j))/(dielecons4solute*sqrt(rij))
                            energy =energy+vdw+ele
                        endif
					enddo
					do j=1, group(k)%cnum2
						rij=(group(ic)%coo2(i,1)-group(k)%coo2(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo2(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo2(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion2(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r2(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(group_para(ic)%dielecons2(i).ge.group_para(k)%dielecons2(j)) then
                                dielecons4solute=group_para(ic)%dielecons2(i)
                            else
                                dielecons4solute=group_para(k)%dielecons2(j)
                            endif
                            ele=(group_para(ic)%charge2(i)*group_para(k)%charge2(j))/(dielecons4solute*sqrt(rij))
                            energy =energy +vdw+ele
                        endif
					enddo
					do j=1, group(k)%cnum3
						rij=(group(ic)%coo2(i,1)-group(k)%coo3(j,1))**2+(group(ic)%coo2(i,2)-group(k)%coo3(j,2))**2+ &
							(group(ic)%coo2(i,3)-group(k)%coo3(j,3))**2
                        if (rij.lt.neighbor_cutoff2) then
                            epsion=sqrt(group_para(ic)%epsion2(i)*group_para(k)%epsion3(j))
                            r0=group_para(ic)%r2(i)+group_para(k)%r3(j)
                            acoeff=epsion*(r0**12)
                            bcoeff=epsion*2*(r0**6)
                            vdw=acoeff/(rij**6)-bcoeff/(rij**3)
                            if(group_para(ic)%dielecons2(i).ge.group_para(k)%dielecons3(j)) then
                                dielecons4solute=group_para(ic)%dielecons2(i)
                            else
                                dielecons4solute=group_para(k)%dielecons3(j)
                            endif
                            ele=(group_para(ic)%charge2(i)*group_para(k)%charge3(j))/(dielecons4solute*sqrt(rij))
                            energy =energy+vdw+ele
                        endif
					enddo
				endif
			endif
		enddo
10		continue
	enddo
    !$OMP END DO 

    !$OMP END PARALLEL

	energy=energy+dihedral_weighting_factor*dihedral_energy

	return
	end  subroutine sidechain_energy

	subroutine binding4sidechainoptimization(group, group_para, ic, numex, inb, numex4, inb4, dihedral_num, dihedral, binding_energy)
	implicit none
	integer*8						:: ligansta, liganend
	integer*8						:: recepsta, recepend
	integer*8						:: natom, atomid(atom_num)
	integer*8    					:: numex(atom_num), inb(atom_num,20), numex4(atom_num), inb4(atom_num,60)
	integer*8						:: ic, dihedral_atom, dihedral_num, ip, jp, kp, lp
	double precision							:: mdcrd(atom_num,3), charge(atom_num), epsion(atom_num)
	double precision							:: r(atom_num), rborn(atom_num), fs(atom_num), dielecons(atom_num)
	character*4						:: lbres(atom_num)

	integer*8						:: i, j, k, flag1, flag2
	double precision							:: rij, vdw, ele, sgb
	double precision							:: binding_vdw, binding_ele, binding_sgb, binding_snp, binding_energy
	double precision							:: p1(3), p2(3), p3(3), p4(3), angle, dihedral_energy, potential
	double precision							:: ipres(gnum), totdevdw(4), totdeele(4), totdesgb(4), totdesnp(4)
	type(groupdetails)				:: group(gnum)
	type(energyparameters)			:: group_para(gnum)
	type(dihedralparameters)		:: dihedral
	double precision, dimension(:), allocatable :: alpha_all, alpha_alone

	natom=0
	flag1=0
	flag2=0
	do i=1,gnum
		ipres(i)=natom
		if(i==ic) dihedral_atom=natom

		do j=1, group(i)%cnum1  !!this loop assigns the various paramters into new labels and tells which is receptor which is ligand
			natom=natom+1
			charge(natom)=group_para(i)%charge1(j)
			epsion(natom)=group_para(i)%epsion1(j)
			r(natom)=group_para(i)%r1(j)
			rborn(natom)=group_para(i)%rborn1(j)
			fs(natom)=group_para(i)%fs1(j)
			dielecons(natom)=group_para(i)%dielecons1(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid1(j)
			mdcrd(natom,1)=group(i)%coo1(j,1)
			mdcrd(natom,2)=group(i)%coo1(j,2)
			mdcrd(natom,3)=group(i)%coo1(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo

		do j=1, group(i)%cnum2
			natom=natom+1
			charge(natom)=group_para(i)%charge2(j)
			epsion(natom)=group_para(i)%epsion2(j)
			r(natom)=group_para(i)%r2(j)
			rborn(natom)=group_para(i)%rborn2(j)
			fs(natom)=group_para(i)%fs2(j)
			dielecons(natom)=group_para(i)%dielecons2(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid2(j)
			mdcrd(natom,1)=group(i)%coo2(j,1)
			mdcrd(natom,2)=group(i)%coo2(j,2)
			mdcrd(natom,3)=group(i)%coo2(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo

		do j=1, group(i)%cnum3
			natom=natom+1
			charge(natom)=group_para(i)%charge3(j)
			epsion(natom)=group_para(i)%epsion3(j)
			r(natom)=group_para(i)%r3(j)
			rborn(natom)=group_para(i)%rborn3(j)
			fs(natom)=group_para(i)%fs3(j)
			dielecons(natom)=group_para(i)%dielecons3(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid3(j)
			mdcrd(natom,1)=group(i)%coo3(j,1)
			mdcrd(natom,2)=group(i)%coo3(j,2)
			mdcrd(natom,3)=group(i)%coo3(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo
	enddo
	recepend=natom	

	allocate (alpha_all(natom))
	!$OMP PARALLEL DO
	do i=1, natom
		call effgbradi(one,natom,i,rborn,fs,mdcrd,alpha_all(i))
	enddo
	!$OMP END PARALLEL DO

	!get effective born radii for peptide (don't need to do receptor since we should have already calculated it by this point)
	allocate (alpha_alone(liganend))
	!$OMP PARALLEL DO
	do i=1, liganend
		call effgbradi(ligansta,liganend,i,rborn,fs,mdcrd,alpha_alone(i))
	enddo
	!$OMP END PARALLEL DO

	totdevdw=0.0
	totdeele=0.0
	totdesgb=0.0
	totdesnp=0.0

	!$OMP PARALLEL PRIVATE(i, j, k, vdw, ele, sgb, rij, flag1)

	!$OMP DO REDUCTION (+ : totdevdw, totdeele, totdesgb, energy_per_res, rec_sgb)
	do i=1, liganend
		do j=i, natom
			vdw=0.0
			ele=0.0
			sgb=0.0
			rij=0.0
			if(j.ne.i) then
				rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
				if (rij.lt.neighbor_cutoff2) then
					do k=1,numex(atomid(i))
						if(atomid(j).eq.inb(atomid(i),k)) goto 40
					enddo
					flag1=0
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag1=1
							goto 45
						endif
					enddo
	45				continue
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					if(flag1==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif
				endif
			endif
40			continue
			call sgbcontri(i,j,rij,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
			if(j.eq.i) sgb=sgb/2.0

			totdevdw(1)=totdevdw(1)+vdw
			totdeele(1)=totdeele(1)+ele
			totdesgb(1)=totdesgb(1)+sgb

			if (j.lt.recepsta) then
				totdevdw(3)=totdevdw(3)+vdw
				totdeele(3)=totdeele(3)+ele
				call sgbcontri(i,j,rij,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
				if(j.eq.i) sgb=sgb/2.0
				totdesgb(3)=totdesgb(3)+sgb
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	deallocate (alpha_all)
	deallocate (alpha_alone)

	totdevdw(4)=totdevdw(1)-totdevdw(3)
	totdeele(4)=totdeele(1)-totdeele(3)
	totdesgb(4)=totdesgb(1)-rec_sgb-totdesgb(3)
	!totdesnp(4)=totdesnp(1)-rec_snp-totdesnp(3)

	binding_vdw=totdevdw(4)
	binding_ele=totdeele(4)
	binding_sgb=totdesgb(4)
	binding_snp=totdesnp(4)+(totdevdw(3)+totdeele(3)+totdesgb(3))*weighting_factor

	dihedral_energy=0.0
	!$OMP PARALLEL PRIVATE(ip, jp, kp, lp, p1,p2,p3,p4,angle,i,j, potential)
	!$OMP DO REDUCTION (+ : dihedral_energy) 
	do i=1, dihedral_num
		ip=dihedral%iph(i)+dihedral_atom; jp=dihedral%jph(i)+dihedral_atom
		kp=dihedral%kph(i)+dihedral_atom; lp=dihedral%lph(i)+dihedral_atom
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		do j=1, dihedral%multiply(i)
			potential=dihedral%pk(i,j)*(1+cosd(dihedral%pn(i,j)*angle-dihedral%phase(i,j)))
			dihedral_energy=dihedral_energy+potential
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	binding_energy=binding_vdw+binding_ele+binding_sgb+binding_snp+dihedral_weighting_factor*dihedral_energy

	return
	end subroutine binding4sidechainoptimization

    subroutine bindingenergy_res_res(i_res, i_atom, j_res, j_atom, group, alpha_alone, alpha_all, charge, dielecons, &
                                    totdesgb, totdevdw, totdeele, mdcrd, numex, numex4, inb, inb4, epsion, r)
    implicit none
    integer*8       :: i_res, i_atom, j_res, j_atom, i, j, k, flag1, flag2, atomid(atom_num)
    type(groupdetails)				:: group(gnum)
    double precision, dimension(:) :: alpha_all, alpha_alone
    double precision              :: vdw,ele,sgb, mdcrd(atom_num,3), rij, epsion(atom_num), r(atom_num)
    double precision              :: dielecons(atom_num), charge(atom_num), totdesgb(4), totdevdw(4), totdeele(4)
    integer*8    					:: numex(atom_num), inb(atom_num,20), numex4(atom_num), inb4(atom_num,60)


	j_atom = i_atom
	do i=1,group(i)%cnum1
		do j=1,group(j)%cnum1
			vdw=0.0
			ele=0.0
			sgb=0.0
			if (i_res.eq.j_res.and.i.eq.j) then
				call sgbcontri(i_atom,j_atom,zero_real,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
				totdesgb(3) = totdesgb(3) + sgb/2.0
				energy_per_res(i_res, 3) = energy_per_res(i_res, 3) - sgb/2.0

				call sgbcontri(i_atom,j_atom,zero_real,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
				totdesgb(1) = totdesgb(1) + sgb/2.0
				energy_per_res(i_res, 3) = energy_per_res(i_res, 3) + sgb/2.0
			else
				rij=(mdcrd(i_atom,1)-mdcrd(j_atom,1))**2+(mdcrd(i_atom,2)-mdcrd(j_atom,2))**2+(mdcrd(i_atom,3)-mdcrd(j_atom,3))**2
				!get energy for isolated ligand
				call sgbcontri(i,j,rij,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
				totdesgb(3)=totdesgb(3)+sgb
				!subtract the sgb values here since binding energy subtracts energy of isolated ligand
				energy_per_res(i_res, 3) = energy_per_res(i_res, 3) - sgb/2.0
				energy_per_res(j_res, 3) = energy_per_res(j_res, 3) - sgb/2.0

				!get energy in bound state (both ligand-receptor and ligand-ligand)
				call sgbcontri(i,j,rij,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
				totdesgb(1)=totdesgb(1)+sgb
				energy_per_res(i_res, 3) = energy_per_res(i_res, 3) + sgb/2.0
				!only store energy if it is a ligand residue
				energy_per_res(j_res, 3) = energy_per_res(j_res, 3) + sgb/2.0

				!now get VDW and ELE energy
				!check if atom pair are in excluded list
				flag1 = 0
				do k=1,numex(atomid(i))
					if(atomid(j).eq.inb(atomid(i),k)) then
						flag1 = 1
						exit
					endif
				enddo

				if (flag1.eq.0) then
					flag2=0
					!check if atom in 1-4 connection
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag1=1
							exit
						endif
					enddo

					!calculate VDW and ELE energy
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					!scale if 1-4 connection
					if(flag1==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif

					!store energy i (only for ligand-receptor interactions, since ligand-ligand and receptor-receptor intertactions will cancel)
					totdevdw(3)=totdevdw(3)+vdw
					totdeele(3)=totdeele(3)+ele
					totdevdw(1)=totdevdw(1)+vdw
					totdeele(1)=totdeele(1)+ele
				endif
			endif
			j_atom = j_atom + 1
		enddo

		do j=1,group(j)%cnum2

		enddo

		do j=1,group(j)%cnum3

		enddo
		i_atom = i_atom + 1
	enddo

    end subroutine bindingenergy_res_res

	subroutine bindingenergy(group, group_para, numex, inb, numex4, inb4, binding_energy, binding_vdw, & 
							binding_ele, binding_sgb, binding_snp)
	implicit none
	integer*8						:: ligansta, liganend
	integer*8						:: recepsta, recepend
	integer*8						:: natom, atomid(atom_num)
	integer*8    					:: numex(atom_num), inb(atom_num,20), numex4(atom_num), inb4(atom_num,60)
	double precision							:: mdcrd(atom_num,3), charge(atom_num), epsion(atom_num)
	double precision							:: r(atom_num), rborn(atom_num), fs(atom_num), dielecons(atom_num)
	character*4						:: lbres(atom_num)
	integer*8						:: i, j, k, flag1, flag2
	double precision							:: rij, vdw, ele, sgb !volume, surfarea
	double precision							:: binding_vdw, binding_ele, binding_sgb, binding_snp, binding_energy
	double precision							:: ipres(gnum), totdevdw(4), totdeele(4), totdesgb(4), totdesnp(4)
	type(groupdetails)				:: group(gnum)
	type(energyparameters)			:: group_para(gnum)
	double precision, dimension(:), allocatable :: alpha_all, alpha_alone
	integer*8						:: residue_atom_limits(sitenum+1), cur_res_i, cur_res_j, per_res_i_flag, per_res_j_flag
	double precision							:: foo

	natom=0
	flag1=0
	flag2=0

	do i=1,gnum
		ipres(i)=natom

		do j=1, group(i)%cnum1  !!this loop assigns the various paramters into new labels and tells which is receptor which is ligand
			natom=natom+1
			charge(natom)=group_para(i)%charge1(j)
			epsion(natom)=group_para(i)%epsion1(j)
			r(natom)=group_para(i)%r1(j)
			rborn(natom)=group_para(i)%rborn1(j)
			fs(natom)=group_para(i)%fs1(j)
			dielecons(natom)=group_para(i)%dielecons1(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid1(j)
			mdcrd(natom,1)=group(i)%coo1(j,1)
			mdcrd(natom,2)=group(i)%coo1(j,2)
			mdcrd(natom,3)=group(i)%coo1(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo

		do j=1, group(i)%cnum2
			natom=natom+1
			charge(natom)=group_para(i)%charge2(j)
			epsion(natom)=group_para(i)%epsion2(j)
			r(natom)=group_para(i)%r2(j)
			rborn(natom)=group_para(i)%rborn2(j)
			fs(natom)=group_para(i)%fs2(j)
			dielecons(natom)=group_para(i)%dielecons2(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid2(j)
			mdcrd(natom,1)=group(i)%coo2(j,1)
			mdcrd(natom,2)=group(i)%coo2(j,2)
			mdcrd(natom,3)=group(i)%coo2(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo

		do j=1, group(i)%cnum3
			natom=natom+1
			charge(natom)=group_para(i)%charge3(j)
			epsion(natom)=group_para(i)%epsion3(j)
			r(natom)=group_para(i)%r3(j)
			rborn(natom)=group_para(i)%rborn3(j)
			fs(natom)=group_para(i)%fs3(j)
			dielecons(natom)=group_para(i)%dielecons3(j)
			atomid(natom)=ipres(i)+group_para(i)%atomid3(j)
			mdcrd(natom,1)=group(i)%coo3(j,1)
			mdcrd(natom,2)=group(i)%coo3(j,2)
			mdcrd(natom,3)=group(i)%coo3(j,3)
			lbres(natom)=group(i)%gtype
			if(i.le.sitenum) then
				if(flag1==0) then
					ligansta=natom
					flag1=1
				endif
			else
				if(flag2==0) then
					liganend=natom-1
					recepsta=natom
					flag2=1
				endif
			endif
		enddo
		if (i.le.sitenum) then
			residue_atom_limits(i) = natom
		endif
	enddo

	recepend=natom
	residue_atom_limits(sitenum+1) = natom !not actually necessary to do this, but just doing it for sake of niceness

	!get effective born radii for complex
	allocate (alpha_all(natom))
	!$OMP PARALLEL DO
	do i=1, natom
		call effgbradi(one,natom,i,rborn,fs,mdcrd,alpha_all(i))
	enddo
	!$OMP END PARALLEL DO

	!get effective born radii for peptide alone 
	allocate (alpha_alone(natom))
	!$OMP PARALLEL DO
	do i=1, liganend
		call effgbradi(ligansta,liganend,i,rborn,fs,mdcrd,alpha_alone(i))
	enddo
	!$OMP END PARALLEL DO

	!get effective born radii for receptor alone 
	!$OMP PARALLEL DO
	do i=recepsta, recepend
		call effgbradi(recepsta,recepend,i,rborn,fs,mdcrd,alpha_alone(i))
	enddo
	!$OMP END PARALLEL DO

	totdevdw=0.0
	totdeele=0.0
	totdesgb=0.0
	totdesnp=0.0

	!$OMP PARALLEL PRIVATE(cur_res_i, per_res_i_flag, cur_res_j, per_res_j_flag, i, j, k, vdw, ele, sgb, rij, flag1)
    cur_res_i = 1
	energy_per_res = 0
	per_res_i_flag = 1
    !$OMP DO REDUCTION (+ : totdevdw, totdeele, totdesgb, energy_per_res, rec_sgb)
	do i=1,natom
		if (per_res_i_flag.eq.1.and.i.gt.residue_atom_limits(cur_res_i)) then
			cur_res_i = cur_res_i + 1
			if (cur_res_i.gt.sitenum) then
				per_res_i_flag = 0
			endif
		endif
		per_res_j_flag = per_res_i_flag
		cur_res_j = cur_res_i
		do j=i, natom
			if (per_res_j_flag.eq.1.and.j.gt.residue_atom_limits(cur_res_j)) then
				cur_res_j = cur_res_j + 1
				if (cur_res_j.gt.sitenum) then
					per_res_j_flag = 0
				endif
			endif
			vdw=0.0
			ele=0.0
			sgb=0.0
			if (i.eq.j) then
				if (i.lt.recepsta) then
					call sgbcontri(i,j,zero_real,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
					totdesgb(3) = totdesgb(3) + sgb/2.0
					energy_per_res(cur_res_i, 3) = energy_per_res(cur_res_i, 3) - sgb/2.0
				else if (rec_calculation_flag.eq.1) then
					call sgbcontri(i,j,zero_real,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
					rec_sgb = rec_sgb + sgb/2.0
				endif
				call sgbcontri(i,j,zero_real,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
				totdesgb(1) = totdesgb(1) + sgb/2.0
				if (i.lt.recepsta) then !only true for peptide
					energy_per_res(cur_res_i, 3) = energy_per_res(cur_res_i, 3) + sgb/2.0
				endif
			else
				rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
                if (rij.lt.neighbor_cutoff2) then
                    if(i.lt.recepsta) then
                        if (j.lt.recepsta) then
                            !get energy for isolated ligand
                            call sgbcontri(i,j,rij,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
                            totdesgb(3)=totdesgb(3)+sgb
                            !subtract the sgb values from per-residue decomposition (binding energy subtracts energy of isolated ligand)
                            energy_per_res(cur_res_i, 3) = energy_per_res(cur_res_i, 3) - sgb/2.0
                            energy_per_res(cur_res_j, 3) = energy_per_res(cur_res_j, 3) - sgb/2.0
                        endif
                        !get energy in bound state (both ligand-receptor and ligand-ligand)
                        call sgbcontri(i,j,rij,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
                        totdesgb(1)=totdesgb(1)+sgb
                        energy_per_res(cur_res_i, 3) = energy_per_res(cur_res_i, 3) + sgb/2.0
                        if (j.lt.recepsta) then
                            !only store energy if it is a ligand residue
                            energy_per_res(cur_res_j, 3) = energy_per_res(cur_res_j, 3) + sgb/2.0
                        endif

                        !now get VDW and ELE energy
                        !check if atom pair are in excluded list
                        do k=1,numex(atomid(i))
                            if(atomid(j).eq.inb(atomid(i),k)) goto 40
                        enddo
                        flag1=0

                        !check if atom in 1-4 connection
                        do k=1,numex4(atomid(i))
                            if(atomid(j).eq.inb4(atomid(i),k)) then
                                flag1=1
                                goto 45
                            endif
                        enddo
    45					continue

                        !calculate VDW and ELE energy
                        call vdwcontri(i,j,rij,epsion,r,vdw)
                        call elecontri(i,j,rij,charge,dielecons,ele)
                        !scale if 1-4 connection
                        if(flag1==1) then
                            vdw=vdw/vdw14_coeff
                            ele=ele/ele14_coeff
                        endif

                        !store energy of peptide-peptide interactions in for peptide in unbound state
                        if (j.lt.recepsta) then
                            totdevdw(3)=totdevdw(3)+vdw
                            totdeele(3)=totdeele(3)+ele
                        endif
						!!store energy of peptide-peptide or peptide-receptor interactions in bound state
                        totdevdw(1)=totdevdw(1)+vdw
                        totdeele(1)=totdeele(1)+ele
                        if (j.ge.recepsta) then  !for per-residue energy decomposition, only add peptide-polymer interactions since peptide-peptide VDW and ELE interactions cancel out in bound vs unbound state
                            energy_per_res(cur_res_i, 1) = energy_per_res(cur_res_i, 1) + vdw/2.0
                            energy_per_res(cur_res_i, 2) = energy_per_res(cur_res_i, 2) + ele/2.0
                        endif
                    else 
                        if (rec_calculation_flag.eq.1) then
                            call sgbcontri(i,j,rij,alpha_alone(i),alpha_alone(j),charge,dielecons,sgb)
                            rec_sgb=rec_sgb+sgb
                        endif
                        call sgbcontri(i,j,rij,alpha_all(i),alpha_all(j),charge,dielecons,sgb)
                        totdesgb(1)=totdesgb(1)+sgb
                    endif
                endif
			endif
40			continue
		enddo
	enddo
	!$OMP END DO
    !$OMP END PARALLEL

	rec_calculation_flag = 0

	deallocate(alpha_all)
	deallocate(alpha_alone)

	!get net binding energy
	totdevdw(4)=totdevdw(1)-totdevdw(3)
	totdeele(4)=totdeele(1)-totdeele(3)
	totdesgb(4)=totdesgb(1)-rec_sgb-totdesgb(3)
	!totdesnp(4)=totdesnp(1)-rec_snp-totdesnp(3)

	!net binding energy per res
	do i=1,sitenum
        energy_per_res(i,4) = sum(energy_per_res(i,1:3))
	enddo

	foo = sum(energy_per_res(:,3))

	binding_vdw = totdevdw(4)
	binding_ele = totdeele(4)
	binding_sgb = totdesgb(4)
	binding_snp = totdesnp(4) + (totdevdw(3) + totdeele(3) + totdesgb(3))*weighting_factor
    binding_energy = binding_vdw + binding_ele + binding_sgb + binding_snp

	return
	end subroutine bindingenergy

	subroutine vdwcontri(x,y,rxy,epsion,r,vdw)
	implicit none
	integer*8					:: x, y
	double precision					:: rxy, vdw
	double precision					:: epsion_xy, r_xy
	double precision					:: acoeff, bcoeff
	double precision					:: epsion(atom_num), r(atom_num)

	epsion_xy=sqrt(epsion(x)*epsion(y))
	r_xy=r(x)+r(y)
	acoeff=epsion_xy*(r_xy**12)
	bcoeff=epsion_xy*2*(r_xy**6)
	vdw=acoeff/(rxy**6)-bcoeff/(rxy**3)

	return
	end subroutine vdwcontri

	subroutine	elecontri(x,y,rxy,charge,dielecons,ele)
	implicit none
	integer*8					:: x, y
	double precision					:: rxy,ele
	double precision					:: qx, qy, dielecons4solute
	double precision					:: charge(atom_num),dielecons(atom_num)

	qx=charge(x)
	qy=charge(y)
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	ele=(qx*qy)/(dielecons4solute*sqrt(rxy))

	return
	end subroutine elecontri

	subroutine	sgbcontri(x,y,rxy,alphax,alphay,charge,dielecons,sgb)
	implicit none
	integer*8					:: x, y
	double precision					:: rxy, sgb
	double precision					:: dielecons4water
	double precision					:: fgb, alphax, alphay
	double precision					:: sgb1, sgb2
	double precision					:: qx, qy, dielecons4solute
	double precision					:: charge(atom_num), dielecons(atom_num)

	dielecons4water=80.0
	qx=charge(x)
	qy=charge(y)
	fgb=sqrt(rxy+alphax*alphay*exp(-rxy/(4*alphax*alphay)))
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	sgb1=(1.0/dielecons4solute)-(1.0/dielecons4water)
	sgb2=qx*qy/fgb
	sgb=-sgb1*sgb2

	return
	end subroutine sgbcontri

	subroutine effgbradi(xstart,xend,x,rborn,fs,mdcrd,alphax)
	implicit none
	integer*8					:: x, xend, xstart
	double precision					:: alpha, beta, gamma
	double precision					:: redborn
	double precision					:: psi, integra, alphax
	double precision					:: rborn(atom_num),fs(atom_num),mdcrd(atom_num,3)

	alpha=0.8
	beta=0.0
	gamma=2.91
	redborn=rborn(x)-0.09
	call areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	psi=integra*redborn
	alphax=1.0/(1.0/redborn-tanh(alpha*psi-beta*psi*psi+gamma*psi*psi*psi)/rborn(x))

	return
	end subroutine effgbradi

	subroutine areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	implicit none
	integer*8					:: x, y
	integer*8					:: xstart, xend
	double precision					:: integra, rxy, sum, redborn
	double precision					:: lxy, uxy
	double precision					:: rborn(atom_num), fs(atom_num), mdcrd(atom_num,3)

	integra=0.0
	do y=xstart, xend
		if(y.ne.x) then
			rxy=(mdcrd(x,1)-mdcrd(y,1))**2+(mdcrd(x,2)-mdcrd(y,2))**2+(mdcrd(x,3)-mdcrd(y,3))**2
			rxy=sqrt(rxy)
			redborn=fs(y)*(rborn(y)-0.09)
			call evalualij(x,rxy,redborn,rborn,lxy)
			call evaluauij(x,rxy,redborn,rborn,uxy)
			sum=(1.0/lxy)-(1.0/uxy)+(1.0/(uxy*uxy)-1.0/(lxy*lxy))*rxy/4.0+log(lxy/uxy)/(2.0*rxy)+ &
				(1.0/(lxy*lxy)-1.0/(uxy*uxy))*redborn*redborn/(4*rxy)
			integra=integra+sum
		endif
	enddo
	integra=integra/2.0

	return
	end	subroutine areafract

	subroutine evalualij(x,rxy,redborn,rborn,lxy)
	implicit none
	integer*8					:: x
	double precision					:: rxy, redborn
	double precision					:: lxy
	double precision					:: rborn(atom_num)

	if(rborn(x).le.(rxy-redborn)) then
		lxy=rxy-redborn
	elseif(rborn(x).le.(rxy+redborn)) then
		lxy=rborn(x)-0.09
	else
		lxy=1.0
	endif

	return
	end subroutine evalualij

	subroutine evaluauij(x,rxy,redborn,rborn,uxy)
	implicit none
	integer*8					:: x
	double precision					:: rxy, redborn
	double precision					:: uxy
	double precision					:: rborn(atom_num)

	if(rborn(x).lt.(rxy+redborn)) then
		uxy=rxy+redborn
	else
		uxy=1.0
	endif

	return
	end subroutine evaluauij

	subroutine arvo(start_res,end_res,group,group_para,volume,surfarea)
	implicit none
	integer*8				:: start_res, end_res, ks, kl, ka, ki, i, j, index, np_test
	double precision					:: sa, volume, surfarea
	double precision					:: rborn(atom_num), mdcrd(atom_num,3)
	parameter (ks=5000,kl=500,ka=500,ki=250000)
	integer*8				:: neighbors_number(ks), index_start(ks), neighbors_indices(ki)
	double precision					:: spheres(ks,4), av(2)
	type(groupdetails)		:: group(gnum)
	type(energyparameters)	:: group_para(gnum)

	!grab coordinates and rborn values group and group_para
	index = 1
	do i=start_res,end_res
		do j=1,group(i)%cnum1
			mdcrd(index,1:3) = group(i)%coo1(j,1:3)
			rborn(index) = group_para(i)%rborn1(j)
			index = index + 1
		enddo
		do j=1,group(i)%cnum2
			mdcrd(index,1:3) = group(i)%coo2(j,1:3)
			rborn(index) = group_para(i)%rborn2(j)
			index =  index + 1
		enddo
		do j=1,group(i)%cnum3
			mdcrd(index,1:3) = group(i)%coo3(j,1:3)
			rborn(index) = group_para(i)%rborn3(j)
			index = index + 1
		enddo
	enddo

	if(index.gt.ks) then
		open(10, file="error.txt", access="append")
			write(10,*) "ks is too small to perform ESURF calculation!"
			write(10,*) "Increase in function then recompile!"
		close(10)
		stop
	endif

	do i=1, index
		spheres(i,2)=mdcrd(i,2)
		spheres(i,3)=mdcrd(i,3)
		spheres(i,4)=rborn(i)+1.40
	enddo

    call make_neighbors(one,index,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,index,ki)

	np_test=0
	do while(np_test.eq.0)
		call NPTest(one,index,spheres,neighbors_number,index_start,neighbors_indices,ks,ki,np_test)
		if(np_test.eq.0) then
			sa=0.324d0
			call spheres_rotation(spheres,ks,index,sa)
		endif
	enddo

	volume=0d0
	surfarea=0d0
	do i=1,index
	    call areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ka,ki,av)
	    volume=volume+av(1)
	    surfarea=surfarea+av(2)
	enddo

	return
	end subroutine arvo

	subroutine make_neighbors(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ns,ki)
	implicit none
	integer*8					:: i1, i2, ns, ks, kl, ki, i, j, neighbors
	integer*8				:: neighbors_number(ks), index_start(ks), neighbors_indices(ki), ind(kl)
	double precision					:: spheres(ks,4)

    index_start(i1)=1
    do i=i1,i2
		call NEIGHBOR(i,spheres,ind,ks,kl,ns,neighbors)
		neighbors_number(i)=neighbors
		if(neighbors_number(i).gt.kl) then
			open(10, file="error.txt", access="append")
				write(10,*) "Too small of kl!"
			close(10)
			stop
		endif
		if (neighbors_number(i).le.0) then
			index_start(i+1)=index_start(i)
		else
			index_start(i+1)=index_start(i)+neighbors_number(i)
			do j=1,neighbors_number(i)
				neighbors_indices(index_start(i)+j-1)=ind(j)
			enddo
		endif
    enddo

	return
	end subroutine make_neighbors

	subroutine NEIGHBOR(i,spheres,ind,ks,kl,ns,neighbors)
	implicit none
	integer*8					:: i, k, kl, ks, ns, neighbors
	integer*8				:: neighbors_num
	double precision					:: xi, yi, zi, ri, rk, dd
	integer*8				:: ind(kl)
	double precision					:: spheres(ks,4)

    neighbors_num=0
    xi=spheres(i,1)
	yi=spheres(i,2)
	zi=spheres(i,3)
	ri=spheres(i,4)
    do k=1,ns
		if (k.ne.i) then
			if(dabs(xi-spheres(k,1)).lt.ri+spheres(k,4)) then
				dd=dsqrt((xi-spheres(k,1))**2+(yi-spheres(k,2))**2+(zi-spheres(k,3))**2)
				rk=spheres(k,4)
				if (dd.lt.ri+rk) then
					if (dd+ri.le.rk) then
						neighbors_num=-1
						exit
					elseif (dd+rk.gt.ri) then
						neighbors_num=neighbors_num+1
						ind(neighbors_num)=k
					endif
				endif
			endif
		endif
      enddo
	neighbors=neighbors_num

	return
	end subroutine NEIGHBOR

	subroutine NPTest(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,ki,np_test)
	implicit none
	integer*8			:: i1, i2, ks, ki, i, k, ink, npt, np_test
	double precision			:: eps_north_pole, dmin, d
	integer*8		:: neighbors_number(ks), index_start(ks), neighbors_indices(ki)
	double precision			:: spheres(ks,4)

    eps_north_pole=1d-4
	dmin=10000d0
	do i=i1,i2
		do k=1,neighbors_number(i)
			ink=neighbors_indices(index_start(i)+k-1)
			d=dabs(dsqrt((spheres(i,1)-spheres(ink,1))**2+(spheres(i,2)-spheres(ink,2))**2 &
     		  +(spheres(i,3)+spheres(i,4)-spheres(ink,3))**2)-spheres(ink,4))
			if (d.lt.dmin) then
				dmin=d
			endif
		enddo
	enddo

    if (dmin.lt.eps_north_pole) then
		npt=0
	else
		npt=1
	endif
 	np_test=npt

    return
	end subroutine NPTest

    subroutine spheres_rotation(spheres,ks,ns,sa)
	implicit none
	integer*8					:: ks, ns, i
	double precision					:: sa, ca, x, z
	double precision					:: spheres(ks,4)

	ca=dsqrt(1d0-sa*sa)
	do i=1,ns
 		x=spheres(i,1)
		z=spheres(i,3)
		spheres(i,1)=ca*x-sa*z
		spheres(i,3)=sa*x+ca*z
	enddo

	return
	end subroutine spheres_rotation

	subroutine areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ka,ki,av)
	implicit none
	integer*8					:: i, ks, kl, ka, ki
	integer*8					:: circles_to_arcs, j, npos, narcs
	integer*8				:: nls
	double precision					:: z1, r1
	integer*8				:: neighbors_number(ks), index_start(ks), neighbors_indices(ki), ind(kl)
	double precision					:: spheres(ks,4), av(2), avi(2)
	double precision					:: circles(kl,4), arcs(ka,3), sphere_local(kl,4)

	if (neighbors_number(i).lt.0) then
		av(1)=0d0
		av(2)=0d0
    elseif (neighbors_number(i).eq.0) then
		av(1)=4d0*pi*spheres(i,4)**3/3.d0
		av(2)=4d0*pi*spheres(i,4)**2
	else
		nls=neighbors_number(i)+1
		ind(1)=i
		do j=1,(nls-1)
			ind(j+1)=neighbors_indices(index_start(i)+j-1)
		enddo
		call local_spheres(spheres,ind,sphere_local,nls,ks,kl)
		av(1)=0d0
		av(2)=0d0
		call make_ts_circles(sphere_local,circles,kl,nls)
		call CirclestoArcs(circles,arcs,kl,nls,ka,circles_to_arcs)
		narcs=circles_to_arcs
		npos=0
		do j=1,(nls-1)
			if (circles(j,4).gt.0) then
				npos=npos+1
			endif
		enddo
		z1=sphere_local(1,3)
		r1=sphere_local(1,4)

		if (npos.gt.0) then
			call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
			av(1)=av(1)+avi(1)
			av(2)=av(2)+avi(2)
         else
			call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
			av(1)=av(1)+avi(1)+4d0*pi*sphere_local(1,4)**3/3d0
			av(2)=av(2)+avi(2)+4d0*pi*sphere_local(1,4)**2
         endif
    endif

	return
	end subroutine areavolume

    subroutine local_spheres(spheres,ind,sphere_local,nls,ks,kl)
	implicit none
	integer*8					:: kl, ks, i, j
	integer*8				:: nls
	integer*8				:: ind(kl)
	double precision					:: spheres(ks,4), sphere_local(kl,4)

	do i=1,nls
		do j=1,4
			sphere_local(i,j)=spheres(ind(i),j)
	    enddo
    enddo

	return
	end subroutine local_spheres

	subroutine make_ts_circles(sphere_local,circles,kl,nls)
	implicit none
	integer*8					:: kl, k
	integer*8				:: nls
	double precision					:: r1, dx, dy, a, b, c, d
	double precision					:: circles(kl,4), sphere_local(kl,4)

	r1=sphere_local(1,4)
	do k=1,(nls-1)
		dx=sphere_local(1,1)-sphere_local(k+1,1)
		dy=sphere_local(1,2)-sphere_local(k+1,2)
		a=dx*dx+dy*dy+(sphere_local(1,3)+r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2
		b=8d0*r1*r1*dx
		c=8d0*r1*r1*dy
		d=4d0*r1*r1*(dx*dx+dy*dy+(sphere_local(1,3)-r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2)
		circles(k,1)=-b/(2d0*a)
		circles(k,2)=-c/(2d0*a)
		circles(k,3)=dsqrt((b*b+c*c-4d0*a*d)/(4d0*a*a))
		if (a.gt.0) then
			circles(k,4)=-1
		else
			circles(k,4)=1
		endif
	enddo

	return
	end subroutine make_ts_circles

	subroutine CirclestoArcs(circles,arcs,kl,nls,ka,circles_to_arcs)
	implicit none
	integer*8					:: kl, ka, number_arc, i, j, k
	integer*8					:: nna, new_arcs, circles_to_arcs
	integer*8				:: nls
	double precision					:: circles(kl,4), arcs(ka,3), arcsnew(ka,3)

    number_arc=0
	if (nls.eq.2) then
		number_arc=1
		arcs(1,1)=1
        arcs(1,2)=0d0
        arcs(1,3)=2d0*pi*circles(1,4)
	else
		do i=1,(nls-1)
			call NewArcs(i,circles,arcsnew,kl,ka,nls,new_arcs)
			nna=new_arcs
			if (nna.gt.0) then
				do j=1,nna
					do k=1,3
						arcs(number_arc+j,k)=arcsnew(j,k)
					enddo
				enddo
				number_arc=number_arc+nna
			endif
		enddo
	endif
	circles_to_arcs=number_arc

	return
	end subroutine CirclestoArcs

	subroutine NewArcs(i,circles,arcsnew,kl,ka,nls,new_arcs)
	implicit none
	integer*8					:: kl, ka, i, j, jj
	integer*8					:: new_arcs
	integer*8					:: circle_in_circle, point_in_circle, delete_equal
	integer*8					:: num_arc, num_angle, number_cond, na
	integer*8				:: nls
	double precision					:: ti, si, ri, t, s, r, d
	double precision					:: a1, a2, b1, b2
	double precision					:: circles(kl,4), arcsnew(ka,3), angles(ka)

	num_arc=0
	num_angle=0
	ti=circles(i,1)
	si=circles(i,2)
	ri=circles(i,3)
	do j=1,(nls-1)
		if (j.ne.i) then
			t=circles(j,1)
			s=circles(j,2)
			r=circles(j,3)
			d=dsqrt((ti-t)**2+(si-s)**2)
			if ( (d.lt.r+ri) .AND. (dabs(r-ri).lt.d) ) then
				call circles_intersection(i,j,circles,kl,a1,a2,b1,b2)
				angles(num_angle+1)=a1
                angles(num_angle+2)=a2
                num_angle=num_angle+2
			endif
		endif
	enddo
	if (num_angle.eq.0) then
		number_cond=0
	    do j=1,(nls-1)
			if (j.ne.i) then
				call CircleinCircle(i,j,circles,kl,circle_in_circle)
				number_cond=number_cond+circle_in_circle
			endif
		enddo
		if (number_cond.eq.(nls-2)) then
			num_arc=1
			arcsnew(1,1)=i
			arcsnew(1,2)=0d0
			arcsnew(1,3)=2d0*pi*circles(i,4)
		endif
	else
		if (circles(i,4).gt.0) then
			call mysort(angles,ka,num_angle)
		else
			call mydsort(angles,ka,num_angle)
		endif
		call DeleteEqual(angles,ka,num_angle,delete_equal)
		na=delete_equal
		num_angle=na
		do j=1,(na-1)
			number_cond=0
			do jj=1,(nls-1)
				if (jj.ne.i) then
					t=ti+ri*dcos((angles(j)+angles(j+1))/2d0)
					s=si+ri*dsin((angles(j)+angles(j+1))/2d0)
					call PointinCircle(t,s,jj,circles,kl,point_in_circle)
					number_cond=number_cond+point_in_circle
				endif
			enddo
			if (number_cond.eq.(nls-2)) then
				num_arc=num_arc+1
				arcsnew(num_arc,1)=i
				arcsnew(num_arc,2)=angles(j)
				arcsnew(num_arc,3)=angles(j+1)-angles(j)
			endif
		enddo
		number_cond=0
		do j=1,(nls-1)
			if (j.ne.i) then
				t=ti+ri*dcos((angles(1)+2d0*pi+angles(na))/2d0)
				s=si+ri*dsin((angles(1)+2d0*pi+angles(na))/2d0)
				call PointinCircle(t,s,j,circles,kl,point_in_circle)
				number_cond=number_cond+point_in_circle
			endif
		enddo
		if (number_cond.eq.(nls-2)) then
			num_arc=num_arc+1
			arcsnew(num_arc,1)=i
			arcsnew(num_arc,2)=angles(na)
			arcsnew(num_arc,3)=angles(1)+circles(i,4)*2d0*pi-angles(na)
		endif
	endif
	new_arcs=num_arc

	return
	end subroutine NewArcs

    subroutine circles_intersection(ic1,ic2,circles,kl,a1,a2,b1,b2)
	implicit none
	integer*8					:: ic1, ic2, kl
	double precision					:: a1, a2, b1, b2
	double precision					:: eps_deltat, t1, s1, r1, t2, s2, r2
	double precision					:: A, B, C, D
	double precision					:: circles(kl,4)

	eps_deltat=1d-12
	t1=circles(ic1,1)
	s1=circles(ic1,2)
	r1=circles(ic1,3)
	t2=circles(ic2,1)
	s2=circles(ic2,2)
	r2=circles(ic2,3)
	if(dabs(t2-t1).lt.eps_deltat) then
		B=((r1*r1-r2*r2)/(s2-s1)-(s2-s1))/2d0
		A=dsqrt(r2*r2-B*B)
		if (B.eq.0) then
			b1=0d0
			b2=pi
		elseif (B.gt.0) then
			b1=datan(dabs(B/A))
			b2=pi-b1
		else
			b1=pi+datan(dabs(B/A))
			b2=3d0*pi-b1
		endif
		B=B+s2-s1
		if (B.eq.0) then
			a1=0d0
			a2=pi
		elseif (B.gt.0) then
			a1=datan(dabs(B/A))
			a2=pi-a1
		else
			a1=pi+datan(dabs(B/A))
			a2=3d0*pi-a1
		endif
	else
		C=((r1*r1-r2*r2-(s2-s1)**2)/(t2-t1)-(t2-t1))/2d0
		D=(s1-s2)/(t2-t1)
		B=(-C*D+dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
		A=C+D*B
		if (A.eq.0) then
			if (B.gt.0) then
				b1=pi/2d0
			else
				b1=-pi/2d0
			endif
		elseif (A.gt.0) then
			b1=datan(B/A)
		else
			b1=pi+datan(B/A)
		endif
		B=B+s2-s1
		A=A+t2-t1
		if (A.eq.0) then
			if (B.gt.0) then
				a1=pi/2d0
			else
				a1=-pi/2d0
			endif
		elseif (A.gt.0) then
			a1=datan(B/A)
		else
			a1=pi+datan(B/A)
		endif
		B=(-C*D-dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
		A=C+D*B
		if (A.eq.0) then
			if (B.gt.0) then
				b2=pi/2d0
			else
				b2=-pi/2d0
			endif
		elseif (A.gt.0) then
			b2=datan(B/A)
		else
			b2=pi+datan(B/A)
		endif
		B=B+s2-s1
		A=A+t2-t1
		if (A.eq.0) then
			if (B.gt.0) then
				a2=pi/2d0
			else
				a2=-pi/2d0
			endif
		elseif (A.gt.0) then
			a2=datan(B/A)
		else
			a2=pi+datan(B/A)
		endif
	endif
	if (a1.lt.0) a1=a1+2d0*pi
	if (a2.lt.0) a2=a2+2d0*pi
	if (b1.lt.0) b1=b1+2d0*pi
	if (b2.lt.0) b2=b2+2d0*pi

    return
	end subroutine circles_intersection

	subroutine CircleinCircle(i,k,circles,kl, circle_in_circle)
	implicit none
	integer*8			::i, k, kl
	integer*8			::circle_in_circle
	double precision			::d
	double precision			::circles(kl,4)

	d=dsqrt((circles(i,1)+circles(i,3)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
	if (d.lt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			circle_in_circle=1
		else
			circle_in_circle=0
		endif
	elseif (d.gt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			circle_in_circle=0
		else
			circle_in_circle=1
		endif
	else
		d=dsqrt((circles(i,1)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
		if (d.lt.circles(k,3)) then
			if (circles(k,4).gt.0) then
				circle_in_circle=1
			else
				circle_in_circle=0
			endif
		else
			if (circles(k,4).gt.0) then
				circle_in_circle=0
			else
				circle_in_circle=1
			endif
		endif
	endif

	return
	end	subroutine CircleinCircle

	subroutine PointinCircle(t,s,k,circles,kl,point_in_circle)
	implicit none
	integer*8					:: k, kl
	integer*8					:: point_in_circle
	double precision					:: t, s, d
	double precision					:: circles(kl,4)

	d=dsqrt((t-circles(k,1))**2+(s-circles(k,2))**2)
	if (d.lt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			point_in_circle=1
		else
			point_in_circle=0
		endif
	else
		if (circles(k,4).gt.0) then
			point_in_circle=0
		else
			point_in_circle=1
		endif
	endif

	return
	end	subroutine PointinCircle

	subroutine mysort(angles,ka,num_angle)
	implicit none
	integer*8					:: ka, num_angle, i, ii, j
	double precision					:: amax
	double precision					:: angles(ka)

	do i=1,(num_angle-1)
		ii=i
		amax=angles(i)
		do j=i+1,num_angle
			if (amax.gt.angles(j)) then
				ii=j
				amax=angles(j)
			endif
 		enddo
		if (ii.ne.i) then
			angles(ii)=angles(i)
			angles(i)=amax
		endif
	enddo

	return
	end	subroutine mysort

	subroutine mydsort(angles,ka,num_angle)
	implicit none
	integer*8					:: ka, num_angle, i, ii, j
	double precision					:: amin
	double precision					:: angles(ka)

	do i=1,(num_angle-1)
		ii=i
		amin=angles(i)
		do j=i+1,num_angle
			if (amin.lt.angles(j)) then
				ii=j
				amin=angles(j)
			endif
 		enddo
		if (ii.ne.i) then
			angles(ii)=angles(i)
			angles(i)=amin
		endif
	enddo

	return
	end subroutine mydsort

	subroutine DeleteEqual(angles,ka,num_angle,delete_equal)
	implicit none
	integer*8					:: ka, num_angle, m, i
	integer*8					:: delete_equal
	double precision					:: eps_angle, angle
	double precision					:: angles(ka), anglesnew(ka)

	eps_angle=1d-12
   	m=1
	angle=angles(1)
	anglesnew(1)=angle
	do i=2,num_angle
		if (dabs(angles(i)-angle).gt.eps_angle) then
			angle=angles(i)
			m=m+1
			anglesnew(m)=angle
		endif
	enddo
	delete_equal=m
	do i=1,m
		angles(i)=anglesnew(i)
	enddo

	return
	end	subroutine DeleteEqual

	subroutine avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
	implicit none
	integer*8					:: kl, ka, narcs, k
	double precision					:: z1, r1, eps_two_pi
	double precision					:: t, s, r, A, B, C, SS, rr
	double precision					:: vIone, vItwo, vIthree
	double precision					:: vJone, vJtwo, vJthree
	double precision					:: delta_vint, delta_aint
	double precision					:: be, al, sb, cb, sa, ca
	double precision					:: a1, a2, b1, b2
	double precision					:: circles(kl,4), arcs(ka,3), avi(2)

	eps_two_pi=1d-12
	avi(1)=0d0
	avi(2)=0d0
	do k=1,narcs
	    t=circles(arcs(k,1),1)
	    s=circles(arcs(k,1),2)
	    r=circles(arcs(k,1),3)
	    A=(4d0*r1*r1+t*t+s*s+r*r)/2d0
	    B=t*r
		C=s*r
		SS=dsqrt(A*A-B*B-C*C)
		rr=r*r-A
		if (dabs(dabs(arcs(k,3))-2d0*pi).lt.eps_two_pi) then
			vIone=2d0*pi/SS
			vItwo=2d0*pi*A/(SS**3)
			vIthree=pi*(2d0*A*A+B*B+C*C)/(SS**5)
			vJone=pi+rr/2d0*vIone
			vJtwo=(vIone+rr*vItwo)/4d0
			vJthree=(vItwo+rr*vIthree)/8d0
			delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
	        delta_aint=2d0*vJone*r1**2
			if (arcs(k,3).lt.0) then
			   delta_vint=-delta_vint
			   delta_aint=-delta_aint
              endif
			avi(1)=avi(1)+delta_vint
			avi(2)=avi(2)+delta_aint
		else
			if (arcs(k,3).lt.0) then
				al=arcs(k,2)+arcs(k,3)
				be=arcs(k,2)
			else
				be=arcs(k,2)+arcs(k,3)
				al=arcs(k,2)
			endif
			vIone=2d0*(pi/2d0-datan((A*dcos((be-al)/2d0)+B*dcos((al+be)/2d0)+C*dsin((al+be)/2d0))/ &
    			  (SS*dsin((be-al)/2d0))))/SS
			sb=dsin(be)
			cb=dcos(be)
			sa=dsin(al)
			ca=dcos(al)
			call Fract(A,B,C,sa,ca,one,a1)
			call Fract(A,B,C,sa,ca,one+1,a2)
			call Fract(A,B,C,sb,cb,one,b1)
			call Fract(A,B,C,sb,cb,one+1,b2)
			vItwo=(b1-a1+A*vIone)/(SS*SS)
			vIthree=(b2-a2+(b1-a1)/A+(2d0*A*A+B*B+C*C)*vItwo/A)/(2d0*S*S)
			vJone=((be-al)+rr*vIone)/2d0
			vJtwo=(vIone+rr*vItwo)/4d0
			vJthree=(vItwo+rr*vIthree)/8d0
			delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
	        delta_aint=2d0*vJone*r1**2
			if (arcs(k,3).lt.0) then
			   delta_vint=-delta_vint
			   delta_aint=-delta_aint
            endif
			avi(1)=avi(1)+delta_vint
			avi(2)=avi(2)+delta_aint
		endif
	enddo

	return
	end	subroutine avintegral

	subroutine Fract(A,B,C,sinphi,cosphi,k,fraction)
	implicit none
	integer*8					:: k
	double precision					:: A, B, C, sinphi, cosphi
	double precision					:: fraction

    fraction=(-B*sinphi+C*cosphi)/(A+B*cosphi+C*sinphi)**k

	return
	end	subroutine Fract

end module energy_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module advancedfunction

	use constant
	use randomgenerator
	use input
!	use pdbfile
	use	mathfunction
	use database
	use transplant
	use energy_calculation

	contains
	subroutine replace_extra_tryptophans(group, groupdata_backup, temp_group)
	implicit none
	integer*8 			:: non_trp_count=0, trp_locs(40), non_trp_locs(40) 
	integer*8 			:: i, j, k, num_over, index_trp=0, index_non_trp, flag, flag1, flag2, rotanum, feedback
	type(groupdetails)	:: group(:), temp_group(:), aa_group(40)
	double precision				:: ran
	character*4			:: group_name_1(3), group_name_2(3), aminoacid_name_1
	type(databackup)	:: groupdata_backup(:)
	
	!calculate number of tryptophans
	trp_count = 0
	do i=1, sitenum
		if (group(i)%gtype == "NTRP".or.group(i)%gtype == "TRP".or.group(i)%gtype == "CTRP") then
			trp_locs(trp_count+1) = i
			trp_count = trp_count + 1
		else
			non_trp_locs(non_trp_count+1) = i
			non_trp_count = non_trp_count + 1
		endif
	enddo

	num_over = trp_count - trp_limit
	temp_group = group
	do i=1,num_over
		!pick a tryptophan randomly and replace with a valine (arbitrarily pick valine since it is small and shouldn't have overlaps and is of same type)
		call ran_gen(ran,zero)
		index_trp=int(ran*trp_count-1.0e-3)+1
		if (index_trp.gt.trp_count) index_trp=trp_count
		
		call ran_gen(ran,zero)
		index_non_trp=int(ran*non_trp_count-1.0e-3)+1
		if (index_non_trp.gt. non_trp_count) index_non_trp=non_trp_count

		flag=0
		call groupinfo(group(trp_locs(index_trp))%gtype, group_name_1, flag1)
		call groupinfo("NVAL", group_name_2, flag2)

		aminoacid_name_1=group_name_2(flag1)

		call findrotamer(trp_locs(index_trp), temp_group, aminoacid_name_1, rotanum, aa_group)

		do j=1, rotanum
			call residue_replace(trp_locs(index_trp), group, groupdata_backup, j, aa_group, temp_group)

			call check_transplant(zero, trp_locs(index_trp), zero, temp_group, feedback)

			if(feedback==1) then
				!once we find a rotamer without steric clashes, accept the transplant and move on after updating tryptophan info
				group=temp_group
				non_trp_locs(non_trp_count+1) = trp_locs(index_trp)
				do k=index_trp, trp_count-1
					trp_locs(k) = trp_locs(k+1)
				enddo
				trp_count = trp_count - 1
				non_trp_count = non_trp_count + 1
				goto 10
			endif
		enddo

		!if we reach this part, then we couldn't replace the tryptophan with valine
		!for now, will have program exit since I don't expect this will happen often.
		!if this is a problem, then will add a solution 
		open(20, file="error.txt", access="append")
			write(20,*) "Unable to replace tryptophan with valine! Terminating program."
			close(20)
		stop

		10 continue

	enddo

	return
	end subroutine replace_extra_tryptophans

    function calc_COM(group, flag)
    implicit none
    type(groupdetails)				:: group(gnum)
    integer                       	:: i,j, flag
    double precision                          :: calc_COM(3),  mass, net_mass
    character*4                     :: atom_type

	if (flag.eq.1) then !calculate peptide center of mass 
		calc_COM = 0.0; net_mass = 0.0
		do i=1, sitenum
			do j=1, group(i)%cnum1
				atom_type = group(i)%atype1(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo1(j,:)
				net_mass = net_mass + mass
				if (group(i)%coo1(j,3).lt.pep_min_z) pep_min_z = group(i)%coo1(j,3)

			enddo

			do j=1, group(i)%cnum2
				atom_type = group(i)%atype2(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo2(j,:)
				net_mass = net_mass + mass
				if (group(i)%coo2(j,3).lt.pep_min_z) pep_min_z = group(i)%coo2(j,3)
			enddo

			do j=1, group(i)%cnum3
				atom_type = group(i)%atype3(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo3(j,:)
				net_mass = net_mass + mass
				if (group(i)%coo3(j,3).lt.pep_min_z) pep_min_z = group(i)%coo3(j,3)
			enddo
			calc_COM = calc_COM/net_mass
    	enddo

	elseif (flag.eq.2) then
		calc_COM = 0.0; net_mass = 0.0
		do i=sitenum+1, gnum
			do j=1, group(i)%cnum1
				atom_type = group(i)%atype1(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo1(j,:)
				net_mass = net_mass + mass
			enddo

			do j=1, group(i)%cnum2
				atom_type = group(i)%atype2(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo2(j,:)
				net_mass = net_mass + mass
			enddo

			do j=1, group(i)%cnum3
				atom_type = group(i)%atype3(j)

				if (atom_type(1:1) == 'H') then
					mass = 1.0
				elseif (atom_type(1:1) == 'C') then
					mass = 12.0
				elseif (atom_type(1:1) == 'N') then
					mass = 14.0
				elseif (atom_type(1:1) == 'O') then
					mass = 16.0
				elseif (atom_type(1:1) == 'S') then
					mass = 32.1
				else
					open(10, file='error.txt')
					write(10,*) "Unknown atom type in center of mass calculation! Terminating program"
					close(10)
					stop
				endif

				calc_COM = calc_COM + mass*group(i)%coo3(j,:)
				net_mass = net_mass + mass
			enddo
    	enddo
	else
		open(10,file='error.txt')
		write(10,*) "Flag passed to 'calc_COM' function was ", flag
		write(10,*) "This value is not accepted - there is a bug in the program!!"
		write(10,*) "Terminating program"
		close(10)
		stop
	endif

	calc_COM = calc_COM/net_mass

	end function calc_COM

    subroutine translate_ligand(group, vector)
    type(groupdetails)				:: group(gnum)
    double precision                          :: vector(3)
    integer*8                       :: i,j

    do i=1,sitenum
        do j=1, group(i)%cnum1
            group(i)%coo1(j,:) = group(i)%coo1(j,:) + vector
        enddo

        do j=1, group(i)%cnum2
            group(i)%coo2(j,:) = group(i)%coo2(j,:) + vector
        enddo

        do j=1, group(i)%cnum3
            group(i)%coo3(j,:) = group(i)%coo3(j,:) + vector
        enddo
    enddo

    end subroutine translate_ligand

	subroutine move_pep_over_surface(group)
	implicit none
	real					:: shift(3), z_shift, delta
	integer					:: res, atom
	type(groupdetails)		:: group(gnum)

	shift = 0
	pep_com = calc_COM(group,1)
	pol_com = calc_COM(group,2)
	shift(1:2) = pol_com(1:2) - pep_com(1:2) !get length to shift peptide so its center of mass aligns with polymer's in x/y directions
	z_shift = target_pep_distance + surface_z - pep_com(3)

	pep_min_z = pep_min_z + z_shift
	if (pep_min_z - surface_z .gt. 2) then
		delta = -1*(2 - (pep_min_z - surface_z))
		pep_min_z = pep_min_z + delta
		z_shift = z_shift + delta
	endif
	shift(3) = z_shift
	pep_com = pep_com + shift

	do res=1,sitenum
		do atom=1, group(res)%cnum1
			group(res)%coo1(atom,:) = group(res)%coo1(atom,:) + shift
		enddo

		do atom=1, group(res)%cnum2
			group(res)%coo2(atom,:) = group(res)%coo2(atom,:) + shift
		enddo

		do atom=1, group(res)%cnum3
			group(res)%coo3(atom,:) = group(res)%coo3(atom,:) + shift
		enddo

	enddo

	return
	end subroutine move_pep_over_surface

    subroutine random_ligand_translation(group)
    implicit none
    type(groupdetails)				:: group(gnum)
    double precision                          :: translate_vec(3), val1, val2, val3
	double precision							:: lims(3,2), ran, mag, com(3)
	integer*4						:: i, sel_index, scalar

	!*******NOTE**********!
	!in future, may want to change this function to translation occurs only in one direction at at time, similar to rotation function

    !generate random numbers between -1 and 1
    call ran_gen(val1,zero)
    val1 = val1*2 - 1
    call ran_gen(val2,zero)
    val2 = val2*2 - 1
    call ran_gen(val3,zero)
    val3 = val3*2 - 1

	!determine current center of mass of peptide (current value may not be accurate since center of mass not calculated after backbone or sequence changes)
	com = calc_COM(group,1)

	!determine max translation in each direction allowed while still being within cutoff of starting center of mass position
	!z-direction is special for surface design, so we have extra control over translation distance in z direction

	lims(1,:)= (/com(1) - start_pep_com(1) + max_dev_from_center, start_pep_com(1) - com(1) + max_dev_from_center/)
	lims(2,:)= (/com(2) - start_pep_com(2) + max_dev_from_center, start_pep_com(2) - com(2) + max_dev_from_center/)
	if (translate_z_flag.eq.1) then
		lims(3,:)= (/(com(3)-surface_z) - min_zdist_from_surface, max_zdist_from_surface - (com(3)-surface_z)/)
	else 
		lims(3,:)= (/com(3) - start_pep_com(3) + max_dev_from_center, start_pep_com(3) - com(3) + max_dev_from_center/)
	endif

	translate_vec = 0; scalar=0
	do i =1,3
		sel_index = 0
		!pick directions of translations, only moving in sladirections where minimum translation distance is permissible
		if (lims(i,1).gt.translate_min .and. lims(1,2).gt.translate_min) then
			!either direction works, pick one randomly
			call ran_gen(ran,zero)
			if (ran.lt.0.5) then !go in negative direction
				sel_index = 1
				scalar = -1
			else !go in positive direction
				sel_index = 2
				scalar = 1
			endif
		elseif(lims(i,1).gt.translate_min) then
			sel_index = 1
			scalar = -1
		elseif(lims(i,2).gt.translate_min) then
			sel_index = 2
			scalar = 1
		endif

		!now calculate translation vector (if possible to move in this direction).
		!This translation is done such that a random distance is picked in negative/postive direction (selected above) between
		!translate_min and the minimum of translate_max and limit of peptide displacement from starting center of mass
		if (sel_index.ne.0) then
			mag = min(translate_max, lims(i,sel_index))
			call ran_gen(ran,zero) 
			translate_vec(i) = translate_min + scalar*ran*(mag-translate_min)
		endif
	enddo

    !translate peptide using translation vector calculated in above loop
    call translate_ligand(group, translate_vec)

    end subroutine random_ligand_translation

	subroutine rotate_ligand(group, dir, theta)
	implicit none
	type(groupdetails)				:: group(gnum)
	integer*4						:: dir, i, j
	double precision							:: theta, s, c, m(3,3), min_z, vector(3), com(3)

	s = sin(theta); c = cos(theta)

	!get rotation matrix
	if (dir.eq.0) then
	m = transpose(reshape((/ one_real,   zero_real,     zero_real, &
								zero_real,     c,          s, &
								zero_real,    -s,            c/), shape(m)))
	elseif(dir.eq.1) then
	m = transpose(reshape((/ s      ,   zero_real,      c, &
								zero_real,  one_real,   zero_real, &
								c,        zero_real,     -s/), shape(m)))
	else
	m = transpose(reshape((/ c    ,   s,      zero_real, &
								-s,        c,      zero_real, &
								zero_real,  zero_real,  one_real/), shape(m)))
	endif

	!get center of mass of peptide
	com =  calc_COM(group, 1) 

	!now rotate the ligand
	do i=1, sitenum
		do j=1,group(i)%cnum1
			group(i)%coo1(j,:) = matmul(m,group(i)%coo1(j,:) - com) + com
		enddo
		do j=1,group(i)%cnum2
			group(i)%coo2(j,:) = matmul(m,group(i)%coo2(j,:) - com) + com
		enddo
		do j=1,group(i)%cnum3
			group(i)%coo3(j,:) = matmul(m,group(i)%coo3(j,:) - com) + com
		enddo
	enddo

	!check for steric clashes with polymer surface in z-direction, and move peptide if necessary (only if translate_z_flag is 1)
	if (translate_z_flag.eq.1) then
		min_z = 1000
		do i=1,sitenum
			do j=1,group(i)%cnum1
				min_z = min(group(i)%coo1(j,3),min_z)
			enddo
			do j=1,group(i)%cnum2
				min_z = min(group(i)%coo2(j,3),min_z)
			enddo
			do j=1,group(i)%cnum3
				min_z = min(group(i)%coo3(j,3),min_z)
			enddo
		enddo
		if ((min_z-surface_z).lt.1.55) then !we did find a steric clash, so move the peptide up
			vector(1:2) = 0
			vector(3) = 1.55 + surface_z - min_z
			call translate_ligand(group, vector)
		endif
	endif

	end subroutine rotate_ligand

    subroutine random_ligand_rotation(group)
    implicit none
    type(groupdetails)				:: group(gnum)
    double precision                          :: ran, theta
    integer*4	                    :: dir

    !pick an axis to rotate about: x, y, or z
    call ran_gen(ran,zero)
    if (ran.lt.0.3333) then
        dir = 0 ! x-direction
    elseif(ran.lt.0.66667) then
        dir = 1 !y-direction
    else
        dir = 2 !z-direction
    endif

    !get rotation angle
    call ran_gen(theta,zero)
    theta = rotate_max*(theta*2 - 1)

	!perform rotation
	call rotate_ligand(group, dir, theta)

    end subroutine random_ligand_rotation

	subroutine scmf_substitution(group, groupdata_backup, sub_cycle, temp_group_1)
	implicit none
	integer*8							:: sub_cycle, i, j, k, l, ic, ic1, ic2, flag, flag1, flag2 
	integer*8							:: positive_number, negative_number, aminoacid_number, rotanum, feedback
	integer*8							:: icpointer(6,gnum), pos_type(sitenum), neg_type(sitenum)
	double precision							:: ran2
	character*4						:: aminoacid_name(10), aminoacid_name_1, group_name_1(3), group_name_2(3)
	type(groupdetails)				:: group(gnum), temp_group_1(gnum), aa_group(40)
	type(databackup)				:: groupdata_backup(:)
	integer*8							:: pos_neg_flag !1 if pos > neg, 0 if pos < neg
	integer*8							:: replace_options(6), replace_count, rep_pos_group, rep_neg_group

	cur_hyd=0
	if(ph_value.le.3.9) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. &
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or. &
				group(ic)%gtype=="TRP".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or. &
				group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or. &
				group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP".or. &
				group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP".or. & 
				group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="GLH".or.group(ic)%gtype=="ASH".or. &
				group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or. &
				group(ic)%gtype=="NGLN".or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="NASH".or. &
				group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN" &
				.or.group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CGLH".or.group(ic)%gtype=="CASH") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA".or. &
				group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA" &
				.or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			endif
		enddo

	elseif(ph_value.gt.3.9.and.ph_value.le.4.3) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. & 
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or. & 
				group(ic)%gtype=="TRP".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or. &
				group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or. & 
				group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP".or. &
				group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP".or. & 
				group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. &
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="GLH".or.group(ic)%gtype=="NSER" &
				.or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN" & 
				.or.group(ic)%gtype=="NGLH".or.group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or. &
				group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CGLH") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA".or. &
				group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA" &
				.or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="ASP".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="CASP") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.gt.4.3.and.ph_value.le.6.0) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. &
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or. &
				group(ic)%gtype=="TRP".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or. &
				group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or. &
				group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET" &
				.or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="HIP".or. &
				group(ic)%gtype=="NARG".or.group(ic)%gtype=="NLYS".or.group(ic)%gtype=="NHIP".or. & 
				group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS".or.group(ic)%gtype=="CHIP") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. &
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or. &
				group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="CSER".or. & 
				group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA".or. & 
				group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA" &
				.or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NGLU".or. & 
				group(ic)%gtype=="NASP".or.group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.gt.6.0.and.ph_value.lt.8.3) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE" & 
				.or.group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR" & 
				.or.group(ic)%gtype=="TRP".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL" & 
				.or.group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE" & 
				.or.group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET" & 
				.or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NARG".or. & 
				group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NSER".or. & 
				group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. & 
				group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or. & 
				group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CHIE") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="CYS".or.group(ic)%gtype=="ALA".or. & 
				group(ic)%gtype=="NPRO".or.group(ic)%gtype=="NCYS".or.group(ic)%gtype=="NALA" &
				.or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CCYS".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="NGLU".or. & 
				group(ic)%gtype=="NASP".or.group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.ge.8.3.and.ph_value.lt.10.1) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. & 
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TYR".or. & 
				group(ic)%gtype=="TRP".or.group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or. & 
				group(ic)%gtype=="NILE".or.group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or. & 
				group(ic)%gtype=="NTYR".or.group(ic)%gtype=="NTRP".or.group(ic)%gtype=="CLEU" &
				.or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE".or.group(ic)%gtype=="CMET" & 
				.or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTYR" &
				.or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NARG".or. & 
				group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NSER".or. & 
				group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. & 
				group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or. & 
				group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CHIE") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NPRO".or. & 
				group(ic)%gtype=="NALA".or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or. & 
				group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP".or.group(ic)%gtype=="NCYT" &
			    .or.group(ic)%gtype=="CGLU".or.group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.ge.10.1.and.ph_value.lt.10.5) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. & 
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP".or. & 
				group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or. & 
				group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP".or. & 
				group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="CMET".or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="LYS".or.group(ic)%gtype=="NARG".or. & 
				group(ic)%gtype=="NLYS".or.group(ic)%gtype=="CARG".or.group(ic)%gtype=="CLYS") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="NSER".or. & 
				group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or. & 
				group(ic)%gtype=="NHIE".or.group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or. & 
				group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CHIE") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NPRO".or. & 
				group(ic)%gtype=="NALA".or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or. & 
				group(ic)%gtype=="TYX".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP"  &
				.or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CGLU".or. & 
				group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.ge.10.5.and.ph_value.lt.12.5) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. & 
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP".or. & 
				group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or. & 
				group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP".or. & 
				group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="CMET".or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="ARG".or.group(ic)%gtype=="NARG".or.group(ic)%gtype=="CARG") then
				cur_hyd(3)=cur_hyd(3)+1
				icpointer(3, cur_hyd(3))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="LYN" &
				.or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or.group(ic)%gtype=="NASN".or. & 
				group(ic)%gtype=="NGLN".or.group(ic)%gtype=="NHIE".or.group(ic)%gtype=="NLYN" &
				.or.group(ic)%gtype=="CSER".or.group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or. & 
				group(ic)%gtype=="CGLN".or.group(ic)%gtype=="CHIE".or.group(ic)%gtype=="CLYN") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NPRO".or. & 
				group(ic)%gtype=="NALA".or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or. & 
				group(ic)%gtype=="TYX".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP"  &
				.or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CGLU".or. & 
				group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo

	elseif(ph_value.ge.12.5) then
		do i=1, sitenum
			ic=i
			if(group(ic)%gtype=="GLY".or.group(ic)%gtype=="NGLY".or.group(ic)%gtype=="CGLY") then
				cur_hyd(1)=cur_hyd(1)+1
				icpointer(1, cur_hyd(1))=ic
			elseif(group(ic)%gtype=="LEU".or.group(ic)%gtype=="VAL".or.group(ic)%gtype=="ILE".or. & 
				group(ic)%gtype=="MET".or.group(ic)%gtype=="PHE".or.group(ic)%gtype=="TRP".or. & 
				group(ic)%gtype=="NLEU".or.group(ic)%gtype=="NVAL".or.group(ic)%gtype=="NILE".or. & 
				group(ic)%gtype=="NMET".or.group(ic)%gtype=="NPHE".or.group(ic)%gtype=="NTRP".or. & 
				group(ic)%gtype=="CLEU".or.group(ic)%gtype=="CVAL".or.group(ic)%gtype=="CILE" &
				.or.group(ic)%gtype=="CMET".or.group(ic)%gtype=="CPHE".or.group(ic)%gtype=="CTRP") then
				cur_hyd(2)=cur_hyd(2)+1
				icpointer(2, cur_hyd(2))=ic
			elseif(group(ic)%gtype=="SER".or.group(ic)%gtype=="THR".or.group(ic)%gtype=="ASN".or. & 
				group(ic)%gtype=="GLN".or.group(ic)%gtype=="HIE".or.group(ic)%gtype=="LYN".or. & 
				group(ic)%gtype=="ARN".or.group(ic)%gtype=="NSER".or.group(ic)%gtype=="NTHR".or. & 
				group(ic)%gtype=="NASN".or.group(ic)%gtype=="NGLN".or.group(ic)%gtype=="NHIE".or. & 
				group(ic)%gtype=="NLYN".or.group(ic)%gtype=="NARN".or.group(ic)%gtype=="CSER".or. & 
				group(ic)%gtype=="CTHR".or.group(ic)%gtype=="CASN".or.group(ic)%gtype=="CGLN".or. & 
				group(ic)%gtype=="CHIE".or.group(ic)%gtype=="CLYN".or.group(ic)%gtype=="CARN") then
				cur_hyd(4)=cur_hyd(4)+1
				icpointer(4, cur_hyd(4))=ic
			elseif(group(ic)%gtype=="PRO".or.group(ic)%gtype=="ALA".or.group(ic)%gtype=="NPRO".or. & 
				group(ic)%gtype=="NALA".or.group(ic)%gtype=="CPRO".or.group(ic)%gtype=="CALA") then
				cur_hyd(5)=cur_hyd(5)+1
				icpointer(5, cur_hyd(5))=ic
			elseif(group(ic)%gtype=="GLU".or.group(ic)%gtype=="ASP".or.group(ic)%gtype=="CYT".or. & 
				group(ic)%gtype=="TYX".or.group(ic)%gtype=="NGLU".or.group(ic)%gtype=="NASP"  &
				.or.group(ic)%gtype=="NCYT".or.group(ic)%gtype=="NTYX".or.group(ic)%gtype=="CGLU".or. & 
				group(ic)%gtype=="CASP".or.group(ic)%gtype=="CCYT".or.group(ic)%gtype=="CTYX") then
				cur_hyd(6)=cur_hyd(6)+1
				icpointer(6, cur_hyd(6))=ic
			endif
		enddo
	endif

	!note the indices don't match cause the group order number doesn't match the min/max order
    below_hyd(1) = min_hyd(1) - cur_hyd(1)
    below_hyd(2) = min_hyd(2) - cur_hyd(2)
    below_hyd(3) = min_hyd(3) - cur_hyd(3)
    below_hyd(4) = min_hyd(4) - cur_hyd(4)
    below_hyd(5) = min_hyd(5) - cur_hyd(5)
    below_hyd(6) = min_hyd(6) - cur_hyd(6)

    above_hyd(1) = cur_hyd(1) - max_hyd(1)
    above_hyd(2) = cur_hyd(2) - max_hyd(2)
    above_hyd(3) = cur_hyd(3) - max_hyd(3)
    above_hyd(4) = cur_hyd(4) - max_hyd(4)
    above_hyd(5) = cur_hyd(5) - max_hyd(5)
    above_hyd(6) = cur_hyd(6) - max_hyd(6)

	positive_number=0
	negative_number=0
	do i=1, 6
		if(above_hyd(i)>0) then
            !below loop arbitrarily swaps pointers to amino acids in groups with more amino acid types than desired to randomize order
			do j=1, (cur_hyd(i)-1)
				call ran_gen(ran2,zero)
				ic1=int(ran2*cur_hyd(i)-1.0e-3)+1
				if(ic1.gt.cur_hyd(i)) ic1=cur_hyd(i)
				call ran_gen(ran2,zero)
				ic2=int(ran2*cur_hyd(i)-1.0e-3)+1
				if(ic2.gt.cur_hyd(i)) ic2=cur_hyd(i)

				k=icpointer(i,ic1)
				icpointer(i,ic1)=icpointer(i,ic2)
				icpointer(i,ic2)=k
			enddo
			do j=1, above_hyd(i)
				pos_type(positive_number+j)=i
			enddo
			positive_number=positive_number+above_hyd(i)

		elseif(below_hyd(i)>0) then
			do j=1, below_hyd(i)
				neg_type(negative_number+j)=i
			enddo
			negative_number=negative_number+below_hyd(i)
		endif
	enddo

	do i=1, (positive_number-1)
		call ran_gen(ran2,zero)
		ic1=int(ran2*positive_number-1.0e-3)+1
		if(ic1.gt.positive_number) ic1=positive_number
		call ran_gen(ran2,zero)
		ic2=int(ran2*positive_number-1.0e-3)+1
		if(ic2.gt.positive_number) ic2=positive_number

		k=pos_type(ic1)
		pos_type(ic1)=pos_type(ic2)
		pos_type(ic2)=k
	enddo
	do i=1, (negative_number-1)
		call ran_gen(ran2,zero)
		ic1=int(ran2*negative_number-1.0e-3)+1
		if(ic1.gt.negative_number) ic1=negative_number
		call ran_gen(ran2,zero)
		ic2=int(ran2*negative_number-1.0e-3)+1
		if(ic2.gt.negative_number) ic2=negative_number

		k=neg_type(ic1)
		neg_type(ic1)=neg_type(ic2)
		neg_type(ic2)=k
	enddo

	!temp_group_1=group
	!sub_cycle=min(positive_number, negative_number)
	if (positive_number >= negative_number) then
		sub_cycle = positive_number
		pos_neg_flag = 1
	else		
		sub_cycle = negative_number
		pos_neg_flag = 0
	endif
	do i=1, sub_cycle
		if (pos_neg_flag.eq.1) then
			rep_pos_group = pos_type(i)
			positive_number = positive_number - 1
			if (negative_number.gt.0) then
				rep_neg_group = neg_type(i)
				call scmf_choose_aminoacid(rep_neg_group, aminoacid_number, aminoacid_name)
				negative_number = negative_number - 1
			else
				!find a suitable amino acid type that can be increased while still being in the hydration limits
				replace_count = 0
				do j=1,6
					if (above_hyd(j) < 0) then
						replace_options(replace_count+1) = j
						replace_count = replace_count + 1
					endif
				enddo
				call ran_gen(ran2,zero)
				rep_neg_group=int(ran2*replace_count-1.0e-3)+1
				if (rep_neg_group.gt.replace_count) rep_neg_group=replace_count
				rep_neg_group = replace_options(rep_neg_group)
				call scmf_choose_aminoacid(rep_neg_group, aminoacid_number, aminoacid_name)
			endif
		else
			rep_neg_group = neg_type(i)
			negative_number = negative_number - 1
			if (positive_number.gt.0) then
				rep_pos_group = pos_type(i)
				call scmf_choose_aminoacid(rep_neg_group, aminoacid_number, aminoacid_name)
				positive_number = positive_number - 1
			else
				!find a suitable amino acid type that can be decreased while still being in the hydration limits
				replace_count = 0
				do j=1,6
					if (below_hyd(j) < 0) then
						replace_options(replace_count+1) = j
						replace_count = replace_count + 1
					endif
				enddo
				call ran_gen(ran2,zero)
				rep_pos_group=int(ran2*replace_count-1.0e-3)+1
				if (rep_pos_group.gt.replace_count) rep_pos_group=replace_count
				rep_pos_group = replace_options(rep_pos_group)
				call scmf_choose_aminoacid(rep_neg_group, aminoacid_number, aminoacid_name)
			endif
		endif

		flag=0
		l=1
		do while(l<=cur_hyd(rep_pos_group))
			do j=1, aminoacid_number
                call groupinfo(group(icpointer(rep_pos_group,l))%gtype, group_name_1, flag1)
                call groupinfo(aminoacid_name(j), group_name_2, flag2)

				aminoacid_name_1=group_name_2(flag1)

				call findrotamer(icpointer(rep_pos_group,l), group, aminoacid_name_1, rotanum, aa_group)

				do k=1, rotanum
					call residue_replace(icpointer(rep_pos_group,l), group, groupdata_backup, k, aa_group, temp_group_1)

					call check_transplant(zero, icpointer(rep_pos_group,l), zero, temp_group_1, feedback)

					if(feedback==1) then
						group=temp_group_1
						flag=1
						goto 10
					endif
				enddo
			enddo
			l=l+1
		enddo

10		continue

		if(flag==1) then
			!now fix icpointers so the amino acid types are correct
			icpointer(rep_neg_group,cur_hyd(rep_neg_group)+1) = icpointer(rep_pos_group,l)
			do k=l, cur_hyd(rep_pos_group)
				icpointer(rep_pos_group,k) = icpointer(rep_pos_group,k+1)
			enddo
			cur_hyd(rep_pos_group)=cur_hyd(rep_pos_group) - 1
			cur_hyd(rep_neg_group)=cur_hyd(rep_neg_group) + 1 

			above_hyd(rep_neg_group) = above_hyd(rep_neg_group) + 1
			below_hyd(rep_neg_group) = below_hyd(rep_neg_group) - 1
	
			above_hyd(rep_pos_group) = above_hyd(rep_pos_group) - 1
			below_hyd(rep_pos_group) = below_hyd(rep_pos_group) + 1
		elseif(flag==0) then
			open(20, file="error.txt", access="append")
				write(20,*) "Unable to make scmf_substitution without a steric clash!" 
				write(20,*) "Issue is with residue ", icpointer(rep_pos_group,l), ": ", temp_group_1(icpointer(rep_pos_group,l))
				write(20,*) "Either adjust the hydration properties to permit this amino acid, or provide a different staring structure!"
			close(20)
			stop
		endif
	enddo

	return
	end subroutine scmf_substitution

	subroutine sidechain_optimization(stage, ic, group, group_para, S_numex, S_inb, S_numex4, S_inb4, score)
	implicit none
	integer*8								:: grade, grade_num(6), monitor(6)
	integer*8								:: i, j, k, ic, account_num, flag, stage, trial_count
	integer*8								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer*8								:: dihedral_num
	double precision								:: delta_chi, cos_angle, sin_angle, error, t, Tenergy_min, Tenergy, score, mag
	double precision								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	double precision								:: h2_denominator, h3_denominator
	type(groupdetails)          		:: group(gnum), Tgroup(gnum)
	type(energyparameters)				:: group_para(gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6), class_old(6), class_new(6), class_min(6)
	type(dihedralparameters)			:: dihedral

	double precision, dimension(:), allocatable 	:: energy_forward, energy_backward
	double precision, dimension(:,:), allocatable 	:: gradient_old, gradient_new, Hessian_old, Hessian_new
	double precision, dimension(:,:), allocatable   :: H2, H3, H31
	double precision, dimension(:,:), allocatable	:: d, y, s, Tchi


	call sidechain_category(ic, group, Iclass, grade, grade_num, index, monitor)
	call dihedralangle_reading(group(ic)%gtype, dihedral_num, dihedral)

	allocate(energy_forward(grade)); allocate(energy_backward(grade))
	allocate(gradient_old(grade,1)); allocate(gradient_new(grade,1))
	allocate(Hessian_old(grade,grade)); allocate(Hessian_new(grade,grade))
	allocate(H2(grade,grade)); allocate(H3(grade,grade)); allocate(H31(grade,1))
	allocate(d(grade,1)); allocate(y(grade,1)); allocate(s(grade,1))

	Tgroup=group
	do i=1, Tgroup(ic)%cnum1
		if(Tgroup(ic)%atype1(i)=="CA") then
			CA(1)=Tgroup(ic)%coo1(i,1); CA(2)=Tgroup(ic)%coo1(i,2); CA(3)=Tgroup(ic)%coo1(i,3)
		endif
	enddo
	s=0.0
	class_new=Iclass

30	continue
	if(stage==0) then
		delta_chi=5
	elseif(stage==1) then
		delta_chi=1
	endif
	cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
		endif
        mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
        rotaxis(1)=rotaxis_x/mag
        rotaxis(2)=rotaxis_y/mag
        rotaxis(3)=rotaxis_z/mag

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo

			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember

			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo

		do j=1, Tgroup(ic)%cnum2
			Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
	enddo

	cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
		endif
        mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
        rotaxis(1)=rotaxis_x/mag
        rotaxis(2)=rotaxis_y/mag
        rotaxis(3)=rotaxis_z/mag

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo

			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember

			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo

		do j=1, Tgroup(ic)%cnum2
			Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_backward(i))
	enddo

	do i=1, grade
		gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)
	enddo

	do i=1, grade
		do j=1, grade
			if(i==j) then
				Hessian_new(i,j)=1
			else
				Hessian_new(i,j)=0
			endif
		enddo
	enddo

	account_num=0
	t=0.0
	do while(.true.)
		Hessian_old=Hessian_new
		gradient_old=gradient_new
		d=-matmul(Hessian_old, gradient_old)

		allocate(Tchi(grade,1))
		trial_count=0
		do while(.true.)
			class_old=Iclass
			Tchi=t*d
			flag=0
			if(stage==0) then
				do i=1, grade
					if(Tchi(i,1).gt.5.0) then
						Tchi(i,1)=5.0
						flag=1
					elseif(Tchi(i,1).lt.(-5.0)) then
						Tchi(i,1)=-5.0
						flag=1
					endif
				enddo
			elseif(stage==1) then
				do i=1, grade
					if(Tchi(i,1).gt.1.0) then
						Tchi(i,1)=1.0
						flag=1
					elseif(Tchi(i,1).lt.(-1.0)) then
						Tchi(i,1)=-1.0
						flag=1
					endif
				enddo
			endif
			if(flag==1) account_num=account_num+1

			s=s+Tchi
			do i=1, grade
				cos_angle=cosd(s(i,1)); sin_angle=sind(s(i,1))
				if(i==1) then
					rotaxis_x=class_old(i)%member(monitor(i),1)-CA(1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-CA(2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-CA(3)
				else
					rotaxis_x=class_old(i)%member(monitor(i),1)-class_old(i-1)%member(monitor(i-1),1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-class_old(i-1)%member(monitor(i-1),2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-class_old(i-1)%member(monitor(i-1),3)
				endif
                mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
                rotaxis(1)=rotaxis_x/mag
                rotaxis(2)=rotaxis_y/mag
                rotaxis(3)=rotaxis_z/mag

				call axisrotation(rotaxis, cos_angle, sin_angle, m)

				do j=(i+1), (grade+1)
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=class_old(j)%member(k,1)-class_old(i)%member(monitor(i),1)
						class_old(j)%member(k,2)=class_old(j)%member(k,2)-class_old(i)%member(monitor(i),2)
						class_old(j)%member(k,3)=class_old(j)%member(k,3)-class_old(i)%member(monitor(i),3)
					enddo

					Tmember=matmul(class_old(j)%member, m)
					class_old(j)%member=Tmember

					do k=1, grade_num(j)
						class_old(j)%member(k,1)=anint((class_old(j)%member(k,1)+class_old(i)%member(monitor(i),1))*1000)/1000
						class_old(j)%member(k,2)=anint((class_old(j)%member(k,2)+class_old(i)%member(monitor(i),2))*1000)/1000
						class_old(j)%member(k,3)=anint((class_old(j)%member(k,3)+class_old(i)%member(monitor(i),3))*1000)/1000
					enddo
				enddo
			enddo

			do j=1, Tgroup(ic)%cnum2
				Tgroup(ic)%coo2(j,1)=class_old(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(ic)%coo2(j,2)=class_old(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(ic)%coo2(j,3)=class_old(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, Tenergy)

			if(t==0.0) then
				Tenergy_min=Tenergy
				class_min=class_old
			else
				if(Tenergy.lt.Tenergy_min) then
					Tenergy_min=Tenergy
					class_min=class_old
					trial_count=trial_count+1
				else
					s=s-Tchi
					goto 10
				endif
			endif
			t=delta_chi
		enddo
10		continue
		deallocate(Tchi)
		if(stage==0) then
			if(account_num.gt.20) goto 20
		elseif(stage==1) then
			if(account_num.gt.40) goto 20
		endif

		if(trial_count==0) goto 20
		error=0.0
		do i=1, grade
			error=error+abs(d(i,1))
		enddo
		if(error.lt.0.01) goto 20

		cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)

			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			endif
            mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
            rotaxis(1)=rotaxis_x/mag
            rotaxis(2)=rotaxis_y/mag
            rotaxis(3)=rotaxis_z/mag

			call axisrotation(rotaxis, cos_angle, sin_angle, m)

			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo

				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember

				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo

			do j=1, Tgroup(ic)%cnum2
				Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
		enddo

		cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			endif
            mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
            rotaxis(1)=rotaxis_x/mag
            rotaxis(2)=rotaxis_y/mag
            rotaxis(3)=rotaxis_z/mag

			call axisrotation(rotaxis, cos_angle, sin_angle, m)

			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo

				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember

				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo

			do j=1, Tgroup(ic)%cnum2
				Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, &
								  S_inb4, dihedral_num, dihedral, energy_backward(i))
		enddo

		do i=1, grade
			gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)
		enddo

		y=gradient_new-gradient_old

		H2=matmul(s, transpose(s))
		h2_denominator=0.0
		do i=1, grade
			h2_denominator=h2_denominator+s(i,1)*y(i,1)
		enddo

		H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
		H31=matmul(Hessian_old, y)
		h3_denominator=0.0
		do i=1, grade
			h3_denominator=h3_denominator+y(i,1)*H31(i,1)
		enddo

		Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator

	enddo
20	continue
	class_new=class_min

	if(stage==0.and.Tenergy_min.lt.100.0) then
		stage=1
		goto 30
	endif

	if(stage==0) then
		score=1000.0
	elseif(stage==1) then
		score=Tenergy
	endif

	Iclass=class_min

	deallocate(energy_forward); deallocate(energy_backward)
	deallocate(gradient_old); deallocate(gradient_new)
	deallocate(Hessian_old); deallocate(Hessian_new)
	deallocate(H2); deallocate(H3); deallocate(H31)
	deallocate(d); deallocate(y); deallocate(s)

	do i=1, group(ic)%cnum2
		group(ic)%coo2(i,1)=Iclass(index(i)%class_No)%member(index(i)%member_No,1)
		group(ic)%coo2(i,2)=Iclass(index(i)%class_No)%member(index(i)%member_No,2)
		group(ic)%coo2(i,3)=Iclass(index(i)%class_No)%member(index(i)%member_No,3)
	enddo

	return
	end	subroutine sidechain_optimization

    subroutine numerical_gradient(Tclass, ic, energy_forward, energy_backward, gradient, class_new, grade, grade_num, Tgroup, &
                                 group_para, W_numex, W_inb, W_numex4, W_inb4, dihedral_num, dihedral, monitor, index)
    implicit none
    double precision								:: delta_chi, cos_angle, sin_angle, m(3,3), Tmember(15,3)
    integer*8                           :: i, j, k, grade, monitor(6), grade_num(6), ic, dihedral_num
    type(groupdetails)					:: Tgroup(gnum)
    type(conformer4sidechain)			:: Tclass(6),  class_new(6)
    double precision								:: rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), CA(3), mag
    type(index4sidechain)				:: index(60)
    type(energyparameters)				:: group_para(gnum)
    type(dihedralparameters)			:: dihedral
    integer*8							:: W_numex(atom_num), W_inb(atom_num,20), W_numex4(atom_num), W_inb4(atom_num, 60)
    double precision               	            :: energy_forward(:), energy_backward(:)
    double precision        	                    :: gradient(:,:)

    delta_chi=1

    !first get the forward difference
	cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
		endif
        mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
        rotaxis(1)=rotaxis_x/mag
        rotaxis(2)=rotaxis_y/mag
        rotaxis(3)=rotaxis_z/mag

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo

			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember

			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo

		do j=1, Tgroup(ic)%cnum2
			Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call binding4sidechainoptimization(Tgroup, group_para, ic, W_numex, W_inb, W_numex4, &
		 W_inb4, dihedral_num, dihedral, energy_forward(i))
	enddo

    !now evaluate the backward value
	cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
		endif
        mag = sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
        rotaxis(1)=rotaxis_x/mag
        rotaxis(2)=rotaxis_y/mag
        rotaxis(3)=rotaxis_z/mag

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo

			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember

			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo

		do j=1, Tgroup(ic)%cnum2
			Tgroup(ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call binding4sidechainoptimization(Tgroup, group_para, ic, W_numex, W_inb, W_numex4, &
		 W_inb4, dihedral_num, dihedral, energy_backward(i))
	enddo

    !numerical gradient is estimated by finite difference between forward and backward value 
	do i=1, grade
		gradient(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)
	enddo

    end subroutine numerical_gradient

	subroutine sidechain_optimization4binding(ic, group, group_para, W_numex, W_inb, W_numex4, W_inb4)
	implicit none
	integer*8							:: grade, grade_num(6), monitor(6)
	integer*8							:: i, j, k, ic, account_num, flag, trial_count
	integer*8							:: W_numex(atom_num), W_inb(atom_num,20), W_numex4(atom_num), W_inb4(atom_num, 60)
	integer*8							:: dihedral_num
	double precision								:: delta_chi, cos_angle, sin_angle, error, t, Tenergy_min, Tenergy
	double precision								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	double precision								:: h2_denominator, h3_denominator
	type(groupdetails)					:: group(gnum), Tgroup(gnum)
	type(energyparameters)				:: group_para(gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6), class_old(6), class_new(6), class_min(6)
	type(dihedralparameters)			:: dihedral

	double precision, dimension(:), allocatable 	:: energy_forward, energy_backward
	double precision, dimension(:,:), allocatable :: gradient_old, gradient_new, Hessian_old, Hessian_new
	double precision, dimension(:,:), allocatable :: H2, H3, H31
	double precision, dimension(:,:), allocatable	:: d, y, s, Tchi


	call sidechain_category(ic, group, Iclass, grade, grade_num, index, monitor)
	call dihedralangle_reading(group(ic)%gtype, dihedral_num, dihedral)

	allocate(energy_forward(grade)); allocate(energy_backward(grade))
	allocate(gradient_old(grade,1)); allocate(gradient_new(grade,1))
	allocate(Hessian_old(grade,grade)); allocate(Hessian_new(grade,grade))
	allocate(H2(grade,grade)); allocate(H3(grade,grade)); allocate(H31(grade,1))
	allocate(d(grade,1)); allocate(y(grade,1)); allocate(s(grade,1))

	Tgroup=group
	do i=1, Tgroup(ic)%cnum1
		if(Tgroup(ic)%atype1(i)=="CA") then
			CA(1)=Tgroup(ic)%coo1(i,1); CA(2)=Tgroup(ic)%coo1(i,2); CA(3)=Tgroup(ic)%coo1(i,3)
		endif
	enddo
	s=0.0
	class_new=Iclass

	delta_chi=1

    !calculate numerical gradient
    call numerical_gradient(Tclass, ic, energy_forward, energy_backward, gradient_new, class_new, grade, grade_num, Tgroup, &
                            group_para, W_numex, W_inb, W_numex4, W_inb4, dihedral_num, dihedral, monitor, index)

    !initialize Hessian as identity matrix
	do i=1, grade
		do j=1, grade
			if(i==j) then
				Hessian_new(i,j)=1
			else
				Hessian_new(i,j)=0
			endif
		enddo
	enddo

	account_num=0
	t=0.0
	do while(.true.)
		Hessian_old=Hessian_new
		gradient_old=gradient_new
		d=-matmul(Hessian_old, gradient_old)

		allocate(Tchi(grade,1))
		trial_count=0
		do while(.true.)
			class_old=Iclass
			Tchi=t*d
			flag=0
			do i=1, grade
				if(Tchi(i,1).gt.1.0) then
					Tchi(i,1)=1.0
					flag=1
				elseif(Tchi(i,1).lt.(-1.0)) then
					Tchi(i,1)=-1.0
					flag=1
				endif
			enddo
			if(flag==1) account_num=account_num+1

			s=s+Tchi
			do i=1, grade
				cos_angle=cosd(s(i,1)); sin_angle=sind(s(i,1))
				if(i==1) then
					rotaxis_x=class_old(i)%member(monitor(i),1)-CA(1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-CA(2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-CA(3)
				else
					rotaxis_x=class_old(i)%member(monitor(i),1)-class_old(i-1)%member(monitor(i-1),1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-class_old(i-1)%member(monitor(i-1),2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-class_old(i-1)%member(monitor(i-1),3)
				endif
                rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
                rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
                rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

				call axisrotation(rotaxis, cos_angle, sin_angle, m)

				do j=(i+1), (grade+1)
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=class_old(j)%member(k,1)-class_old(i)%member(monitor(i),1)
						class_old(j)%member(k,2)=class_old(j)%member(k,2)-class_old(i)%member(monitor(i),2)
						class_old(j)%member(k,3)=class_old(j)%member(k,3)-class_old(i)%member(monitor(i),3)
					enddo

					Tmember=matmul(class_old(j)%member, m)
					class_old(j)%member=Tmember

					do k=1, grade_num(j)
						class_old(j)%member(k,1)=anint((class_old(j)%member(k,1)+class_old(i)%member(monitor(i),1))*1000)/1000
						class_old(j)%member(k,2)=anint((class_old(j)%member(k,2)+class_old(i)%member(monitor(i),2))*1000)/1000
						class_old(j)%member(k,3)=anint((class_old(j)%member(k,3)+class_old(i)%member(monitor(i),3))*1000)/1000
					enddo
				enddo
			enddo

			do j=1, Tgroup(ic)%cnum2
				Tgroup(ic)%coo2(j,1)=class_old(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(ic)%coo2(j,2)=class_old(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(ic)%coo2(j,3)=class_old(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call binding4sidechainoptimization(Tgroup, group_para, ic, W_numex, W_inb, W_numex4, W_inb4, dihedral_num, dihedral, Tenergy)

			if(t==0.0) then
				Tenergy_min=Tenergy
				class_min=class_old
			else
				if(Tenergy.lt.Tenergy_min) then
					Tenergy_min=Tenergy
					class_min=class_old
					trial_count=trial_count+1
				else
					s=s-Tchi
					goto 10
				endif
			endif
			t=delta_chi
		enddo
10		continue
		deallocate(Tchi)
		if(account_num.gt.40) goto 20

		if(trial_count==0) goto 20
		error=0.0
		do i=1, grade
			error=error+abs(d(i,1))
		enddo
		if(error.lt.0.01) goto 20

        !recalculate numerical gradient
        call numerical_gradient(Tclass, ic, energy_forward, energy_backward, gradient_new, class_min, grade, grade_num, Tgroup, &
        group_para, W_numex, W_inb, W_numex4, W_inb4, dihedral_num, dihedral, monitor,index)

        !calculate change in gradient
		y=gradient_new-gradient_old

		H2=matmul(s, transpose(s))
		h2_denominator=0.0
		do i=1, grade
			h2_denominator=h2_denominator+s(i,1)*y(i,1)
		enddo

		H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
		H31=matmul(Hessian_old, y)
		h3_denominator=0.0
		do i=1, grade
			h3_denominator=h3_denominator+y(i,1)*H31(i,1)
		enddo

		Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator

	enddo
20	continue
	Iclass=class_min

	deallocate(energy_forward); deallocate(energy_backward)
	deallocate(gradient_old); deallocate(gradient_new)
	deallocate(Hessian_old); deallocate(Hessian_new)
	deallocate(H2); deallocate(H3); deallocate(H31)
	deallocate(d); deallocate(y); deallocate(s)

	do i=1, group(ic)%cnum2
		group(ic)%coo2(i,1)=Iclass(index(i)%class_No)%member(index(i)%member_No,1)
		group(ic)%coo2(i,2)=Iclass(index(i)%class_No)%member(index(i)%member_No,2)
		group(ic)%coo2(i,3)=Iclass(index(i)%class_No)%member(index(i)%member_No,3)
	enddo

	return
	end	subroutine sidechain_optimization4binding

	subroutine backbone_rotation_center(group, groupdata_backup, ran_resi, phi1, psi1, &
										 phi2, psi2, Tgroup, Tgroupdata_backup)
	implicit none
	integer*8					:: ran_resi(3), i, j, k, l, temp_num
	double precision					:: C0(3), N1(3), CA1(3), C1(3), N2(3), CA2(3), C2(3), N3(3), N(3), CA(3), C(3)
	double precision					:: rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3)
	double precision					:: phi_old(2), psi_old(2), phi_new(2), psi_new(2), phi1, psi1, phi2, psi2
	double precision					:: delta_phi, delta_psi, cos_angle, sin_angle
	double precision					:: temp1(20,3), temp2(60,3), temp3(20,3), coo(20,3)
	double precision					:: T_backup_temp(3)
	character*4				:: atype(20)

	type(groupdetails)		:: group(gnum), Tgroup(gnum)
	type(databackup)		:: groupdata_backup(gnum), Tgroupdata_backup(gnum)

	phi_new(1)=phi1; phi_new(2)=phi2
	psi_new(1)=psi1; psi_new(2)=psi2

	do i=1, group(ran_resi(1)-1)%cnum3
		if(group(ran_resi(1)-1)%atype3(i)=="C") then
			C0(1)=group(ran_resi(1)-1)%coo3(i,1)
			C0(2)=group(ran_resi(1)-1)%coo3(i,2)
			C0(3)=group(ran_resi(1)-1)%coo3(i,3)
		endif
	enddo

	do i=1, group(ran_resi(1))%cnum1
		if(group(ran_resi(1))%atype1(i)=="N") then
			N1(1)=group(ran_resi(1))%coo1(i,1)
			N1(2)=group(ran_resi(1))%coo1(i,2)
			N1(3)=group(ran_resi(1))%coo1(i,3)
		elseif(group(ran_resi(1))%atype1(i)=="CA") then
			CA1(1)=group(ran_resi(1))%coo1(i,1)
			CA1(2)=group(ran_resi(1))%coo1(i,2)
			CA1(3)=group(ran_resi(1))%coo1(i,3)
		endif
	enddo
	do i=1, group(ran_resi(1))%cnum3
		if(group(ran_resi(1))%atype3(i)=="C") then
			C1(1)=group(ran_resi(1))%coo3(i,1)
			C1(2)=group(ran_resi(1))%coo3(i,2)
			C1(3)=group(ran_resi(1))%coo3(i,3)
		endif
	enddo

	do i=1, group(ran_resi(2))%cnum1
		if(group(ran_resi(2))%atype1(i)=="N") then
			N2(1)=group(ran_resi(2))%coo1(i,1)
			N2(2)=group(ran_resi(2))%coo1(i,2)
			N2(3)=group(ran_resi(2))%coo1(i,3)
		elseif(group(ran_resi(2))%atype1(i)=="CA") then
			CA2(1)=group(ran_resi(2))%coo1(i,1)
			CA2(2)=group(ran_resi(2))%coo1(i,2)
			CA2(3)=group(ran_resi(2))%coo1(i,3)
		endif
	enddo
	do i=1, group(ran_resi(2))%cnum3
		if(group(ran_resi(2))%atype3(i)=="C") then
			C2(1)=group(ran_resi(2))%coo3(i,1)
			C2(2)=group(ran_resi(2))%coo3(i,2)
			C2(3)=group(ran_resi(2))%coo3(i,3)
		endif
	enddo

	do i=1, group(ran_resi(3))%cnum1
		if(group(ran_resi(3))%atype1(i)=="N") then
			N3(1)=group(ran_resi(3))%coo1(i,1)
			N3(2)=group(ran_resi(3))%coo1(i,2)
			N3(3)=group(ran_resi(3))%coo1(i,3)
		endif
	enddo

	call phipsiomg_angle(C0, N1, CA1, C1, phi_old(1))
	call phipsiomg_angle(N1, CA1, C1, N2, psi_old(1))
	call phipsiomg_angle(C1, N2, CA2, C2, phi_old(2))
	call phipsiomg_angle(N2, CA2, C2, N3, psi_old(2))

	Tgroup=group
	Tgroupdata_backup=groupdata_backup

	do i=1, 2
		delta_phi=-(phi_new(i)-phi_old(i))

		cos_angle=cosd(delta_phi)
		sin_angle=sind(delta_phi)

		do k=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(k)=="N") then
				N(1)=Tgroup(ran_resi(i))%coo1(k,1)
				N(2)=Tgroup(ran_resi(i))%coo1(k,2)
				N(3)=Tgroup(ran_resi(i))%coo1(k,3)
			elseif(Tgroup(ran_resi(i))%atype1(k)=="CA") then
				CA(1)=Tgroup(ran_resi(i))%coo1(k,1)
				CA(2)=Tgroup(ran_resi(i))%coo1(k,2)
				CA(3)=Tgroup(ran_resi(i))%coo1(k,3)
			endif
		enddo

		rotaxis_x=CA(1)-N(1)
		rotaxis_y=CA(2)-N(2)
		rotaxis_z=CA(3)-N(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		temp_num=0
		do j=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(j)=="HA".or.Tgroup(ran_resi(i))%atype1(j)=="HA2".or. &
			Tgroup(ran_resi(i))%atype1(j)=="HA3") then
				temp_num=temp_num+1
				atype(temp_num)=Tgroup(ran_resi(i))%atype1(j)
				coo(temp_num,1)=Tgroup(ran_resi(i))%coo1(j,1)-N(1)
				coo(temp_num,2)=Tgroup(ran_resi(i))%coo1(j,2)-N(2)
				coo(temp_num,3)=Tgroup(ran_resi(i))%coo1(j,3)-N(3)
			endif
		enddo

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=Tgroup(ran_resi(i))%coo2(j,1)-N(1)
			Tgroup(ran_resi(i))%coo2(j,2)=Tgroup(ran_resi(i))%coo2(j,2)-N(2)
			Tgroup(ran_resi(i))%coo2(j,3)=Tgroup(ran_resi(i))%coo2(j,3)-N(3)
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=Tgroup(ran_resi(i))%coo3(j,1)-N(1)
			Tgroup(ran_resi(i))%coo3(j,2)=Tgroup(ran_resi(i))%coo3(j,2)-N(2)
			Tgroup(ran_resi(i))%coo3(j,3)=Tgroup(ran_resi(i))%coo3(j,3)-N(3)
		enddo

		temp1=matmul(coo, m)
		coo=temp1

		temp2=matmul(Tgroup(ran_resi(i))%coo2, m)
		Tgroup(ran_resi(i))%coo2=temp2

		temp3=matmul(Tgroup(ran_resi(i))%coo3, m)
		Tgroup(ran_resi(i))%coo3=temp3

		do l=1, temp_num
			do j=1, Tgroup(ran_resi(i))%cnum1
				if(atype(l)==Tgroup(ran_resi(i))%atype1(j)) then
					Tgroup(ran_resi(i))%coo1(j,1)=anint((coo(l,1)+N(1))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,2)=anint((coo(l,2)+N(2))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,3)=anint((coo(l,3)+N(3))*1000)/1000
				endif
			enddo
		enddo

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=anint((Tgroup(ran_resi(i))%coo2(j,1)+N(1))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,2)=anint((Tgroup(ran_resi(i))%coo2(j,2)+N(2))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,3)=anint((Tgroup(ran_resi(i))%coo2(j,3)+N(3))*1000)/1000
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=anint((Tgroup(ran_resi(i))%coo3(j,1)+N(1))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,2)=anint((Tgroup(ran_resi(i))%coo3(j,2)+N(2))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,3)=anint((Tgroup(ran_resi(i))%coo3(j,3)+N(3))*1000)/1000
		enddo

		do j=ran_resi(i)+1, ran_resi(3)
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-N(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-N(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-N(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-N(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-N(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-N(3)
			enddo
			if(j.ne.ran_resi(3)) then
				do l=1, Tgroup(j)%cnum3
					Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-N(1)
					Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-N(2)
					Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-N(3)
				enddo
			endif
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-N(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-N(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-N(3)

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			if(j.ne.ran_resi(3)) then
				temp3=matmul(Tgroup(j)%coo3, m)
				Tgroup(j)%coo3=temp3
			endif

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+N(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+N(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+N(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+N(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+N(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+N(3))*1000)/1000
			enddo
			if(j.ne.ran_resi(3)) then
				do l=1, Tgroup(j)%cnum3
					Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+N(1))*1000)/1000
					Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+N(2))*1000)/1000
					Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+N(3))*1000)/1000
				enddo
			endif
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+N(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+N(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+N(3))*1000)/1000
		enddo

		delta_psi=-(psi_new(i)-psi_old(i))

		cos_angle=cosd(delta_psi)
		sin_angle=sind(delta_psi)

		do k=1, Tgroup(ran_resi(i))%cnum3
			if(Tgroup(ran_resi(i))%atype3(k)=="C") then
				C(1)=Tgroup(ran_resi(i))%coo3(k,1)
				C(2)=Tgroup(ran_resi(i))%coo3(k,2)
				C(3)=Tgroup(ran_resi(i))%coo3(k,3)
			endif
		enddo

		rotaxis_x=C(1)-CA(1)
		rotaxis_y=C(2)-CA(2)
		rotaxis_z=C(3)-CA(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=Tgroup(ran_resi(i))%coo3(j,1)-CA(1)
			Tgroup(ran_resi(i))%coo3(j,2)=Tgroup(ran_resi(i))%coo3(j,2)-CA(2)
			Tgroup(ran_resi(i))%coo3(j,3)=Tgroup(ran_resi(i))%coo3(j,3)-CA(3)
		enddo

		temp3=matmul(Tgroup(ran_resi(i))%coo3, m)
		Tgroup(ran_resi(i))%coo3=temp3

		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=anint((Tgroup(ran_resi(i))%coo3(j,1)+CA(1))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,2)=anint((Tgroup(ran_resi(i))%coo3(j,2)+CA(2))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,3)=anint((Tgroup(ran_resi(i))%coo3(j,3)+CA(3))*1000)/1000
		enddo

		do j=ran_resi(i)+1, ran_resi(3)
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-CA(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-CA(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-CA(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-CA(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-CA(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-CA(3)
			enddo
			if(j.ne.ran_resi(3)) then
				do l=1, Tgroup(j)%cnum3
					Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-CA(1)
					Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-CA(2)
					Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-CA(3)
				enddo
			endif
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-CA(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-CA(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-CA(3)

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			if(j.ne.ran_resi(3)) then
				temp3=matmul(Tgroup(j)%coo3, m)
				Tgroup(j)%coo3=temp3
			endif

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+CA(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+CA(3))*1000)/1000
			enddo
			if(j.ne.ran_resi(3)) then
				do l=1, Tgroup(j)%cnum3
					Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+CA(1))*1000)/1000
					Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+CA(2))*1000)/1000
					Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+CA(3))*1000)/1000
				enddo
			endif
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+CA(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+CA(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+CA(3))*1000)/1000
		enddo
	enddo

	return
	end subroutine backbone_rotation_center

	subroutine backbone_rotation_Nterm(group, groupdata_backup, ran_resi, delta_psi3, delta_phi3, &
										delta_psi2, delta_phi2, Tgroup, Tgroupdata_backup)
	implicit none
	integer*8				:: ran_resi(3), i, j, k, l, temp_num, ic
	double precision					:: ran2
	double precision					:: CA1(3), C1(3), N(3), CA(3), C(3)
	double precision					:: rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3)
	double precision					:: delta_psi3, delta_phi3, delta_psi2, delta_phi2, delta_psi1
	double precision					:: delta_psi(3), delta_phi(3), cos_angle, sin_angle
	double precision					:: temp1(20,3), temp2(60,3), temp3(20,3), coo(20,3)
	double precision					:: T_backup_temp(3)
	character*4				:: atype(20)
	type(groupdetails)		:: group(gnum), Tgroup(gnum)
	type(databackup)		:: groupdata_backup(gnum), Tgroupdata_backup(gnum)

	delta_psi(3)=-delta_psi3; delta_psi(2)=-delta_psi2
	delta_phi(3)=-delta_phi3; delta_phi(2)=-delta_phi2

	Tgroup=group
	Tgroupdata_backup=groupdata_backup

	do i=3, 2, -1
		cos_angle=cosd(delta_psi(i))
		sin_angle=sind(delta_psi(i))

		do k=1, Tgroup(ran_resi(i))%cnum3
			if(Tgroup(ran_resi(i))%atype3(k)=="C") then
				C(1)=Tgroup(ran_resi(i))%coo3(k,1)
				C(2)=Tgroup(ran_resi(i))%coo3(k,2)
				C(3)=Tgroup(ran_resi(i))%coo3(k,3)
			endif
		enddo
		do k=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(k)=="CA") then
				CA(1)=Tgroup(ran_resi(i))%coo1(k,1)
				CA(2)=Tgroup(ran_resi(i))%coo1(k,2)
				CA(3)=Tgroup(ran_resi(i))%coo1(k,3)
			endif
		enddo

		rotaxis_x=CA(1)-C(1)
		rotaxis_y=CA(2)-C(2)
		rotaxis_z=CA(3)-C(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=Tgroup(ran_resi(i))%coo2(j,1)-C(1)
			Tgroup(ran_resi(i))%coo2(j,2)=Tgroup(ran_resi(i))%coo2(j,2)-C(2)
			Tgroup(ran_resi(i))%coo2(j,3)=Tgroup(ran_resi(i))%coo2(j,3)-C(3)
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum1
			Tgroup(ran_resi(i))%coo1(j,1)=Tgroup(ran_resi(i))%coo1(j,1)-C(1)
			Tgroup(ran_resi(i))%coo1(j,2)=Tgroup(ran_resi(i))%coo1(j,2)-C(2)
			Tgroup(ran_resi(i))%coo1(j,3)=Tgroup(ran_resi(i))%coo1(j,3)-C(3)
		enddo
		Tgroupdata_backup(ran_resi(i))%coo(1)=Tgroupdata_backup(ran_resi(i))%coo(1)-C(1)
		Tgroupdata_backup(ran_resi(i))%coo(2)=Tgroupdata_backup(ran_resi(i))%coo(2)-C(2)
		Tgroupdata_backup(ran_resi(i))%coo(3)=Tgroupdata_backup(ran_resi(i))%coo(3)-C(3)

		temp2=matmul(Tgroup(ran_resi(i))%coo2, m)
		Tgroup(ran_resi(i))%coo2=temp2

		temp1=matmul(Tgroup(ran_resi(i))%coo1, m)
		Tgroup(ran_resi(i))%coo1=temp1

		T_backup_temp=matmul(Tgroupdata_backup(ran_resi(i))%coo, m)
		Tgroupdata_backup(ran_resi(i))%coo=T_backup_temp

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=anint((Tgroup(ran_resi(i))%coo2(j,1)+C(1))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,2)=anint((Tgroup(ran_resi(i))%coo2(j,2)+C(2))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,3)=anint((Tgroup(ran_resi(i))%coo2(j,3)+C(3))*1000)/1000
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum1
			Tgroup(ran_resi(i))%coo1(j,1)=anint((Tgroup(ran_resi(i))%coo1(j,1)+C(1))*1000)/1000
			Tgroup(ran_resi(i))%coo1(j,2)=anint((Tgroup(ran_resi(i))%coo1(j,2)+C(2))*1000)/1000
			Tgroup(ran_resi(i))%coo1(j,3)=anint((Tgroup(ran_resi(i))%coo1(j,3)+C(3))*1000)/1000
		enddo
		Tgroupdata_backup(ran_resi(i))%coo(1)=anint((Tgroupdata_backup(ran_resi(i))%coo(1)+C(1))*1000)/1000
		Tgroupdata_backup(ran_resi(i))%coo(2)=anint((Tgroupdata_backup(ran_resi(i))%coo(2)+C(2))*1000)/1000
		Tgroupdata_backup(ran_resi(i))%coo(3)=anint((Tgroupdata_backup(ran_resi(i))%coo(3)+C(3))*1000)/1000

		do j=ran_resi(i)-1, ran_resi(1), -1
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-C(1)
				Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-C(2)
				Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-C(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-C(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-C(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-C(3)
			enddo
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-C(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-C(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-C(3)
			enddo
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-C(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-C(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-C(3)

			temp3=matmul(Tgroup(j)%coo3, m)
			Tgroup(j)%coo3=temp3

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+C(1))*1000)/1000
				Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+C(2))*1000)/1000
				Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+C(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+C(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+C(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+C(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+C(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+C(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+C(3))*1000)/1000
			enddo
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+C(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+C(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+C(3))*1000)/1000
		enddo

		cos_angle=cosd(delta_phi(i))
		sin_angle=sind(delta_phi(i))

		do k=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(k)=="N") then
				N(1)=Tgroup(ran_resi(i))%coo1(k,1)
				N(2)=Tgroup(ran_resi(i))%coo1(k,2)
				N(3)=Tgroup(ran_resi(i))%coo1(k,3)
			endif
		enddo

		rotaxis_x=N(1)-CA(1)
		rotaxis_y=N(2)-CA(2)
		rotaxis_z=N(3)-CA(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		temp_num=0
		do j=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(j)=="N".or.Tgroup(ran_resi(i))%atype1(j)=="H".or.& 
			Tgroup(ran_resi(i))%atype1(j)=="H1".or. &
			Tgroup(ran_resi(i))%atype1(j)=="H2".or.Tgroup(ran_resi(i))%atype1(j)=="H3") then
				temp_num=temp_num+1
				atype(temp_num)=Tgroup(ran_resi(i))%atype1(j)
				coo(temp_num,1)=Tgroup(ran_resi(i))%coo1(j,1)-CA(1)
				coo(temp_num,2)=Tgroup(ran_resi(i))%coo1(j,2)-CA(2)
				coo(temp_num,3)=Tgroup(ran_resi(i))%coo1(j,3)-CA(3)
			endif
		enddo
		Tgroupdata_backup(ran_resi(i))%coo(1)=Tgroupdata_backup(ran_resi(i))%coo(1)-CA(1)
		Tgroupdata_backup(ran_resi(i))%coo(2)=Tgroupdata_backup(ran_resi(i))%coo(2)-CA(2)
		Tgroupdata_backup(ran_resi(i))%coo(3)=Tgroupdata_backup(ran_resi(i))%coo(3)-CA(3)

		temp1=matmul(coo, m)
		coo=temp1

		T_backup_temp=matmul(Tgroupdata_backup(ran_resi(i))%coo, m)
		Tgroupdata_backup(ran_resi(i))%coo=T_backup_temp

		do l=1, temp_num
			do j=1, Tgroup(ran_resi(i))%cnum1
				if(atype(l)==Tgroup(ran_resi(i))%atype1(j)) then
					Tgroup(ran_resi(i))%coo1(j,1)=anint((coo(l,1)+CA(1))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,2)=anint((coo(l,2)+CA(2))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,3)=anint((coo(l,3)+CA(3))*1000)/1000
				endif
			enddo
		enddo
		Tgroupdata_backup(ran_resi(i))%coo(1)=anint((Tgroupdata_backup(ran_resi(i))%coo(1)+CA(1))*1000)/1000
		Tgroupdata_backup(ran_resi(i))%coo(2)=anint((Tgroupdata_backup(ran_resi(i))%coo(2)+CA(2))*1000)/1000
		Tgroupdata_backup(ran_resi(i))%coo(3)=anint((Tgroupdata_backup(ran_resi(i))%coo(3)+CA(3))*1000)/1000

		do j=ran_resi(i)-1, ran_resi(1), -1
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-CA(1)
				Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-CA(2)
				Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-CA(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-CA(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-CA(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-CA(3)
			enddo
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-CA(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-CA(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-CA(3)
			enddo
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-CA(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-CA(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-CA(3)

			temp3=matmul(Tgroup(j)%coo3, m)
			Tgroup(j)%coo3=temp3

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+CA(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+CA(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+CA(3))*1000)/1000
			enddo
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+CA(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+CA(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+CA(3))*1000)/1000
		enddo
	enddo

	do i=1, Tgroup(ran_resi(1))%cnum3
		if(Tgroup(ran_resi(1))%atype3(i)=="C") then
			C1(1)=Tgroup(ran_resi(1))%coo3(i,1)
			C1(2)=Tgroup(ran_resi(1))%coo3(i,2)
			C1(3)=Tgroup(ran_resi(1))%coo3(i,3)
		endif
	enddo
	do i=1, Tgroup(ran_resi(1))%cnum1
		if(Tgroup(ran_resi(1))%atype1(i)=="CA") then
			CA1(1)=Tgroup(ran_resi(1))%coo1(i,1)
			CA1(2)=Tgroup(ran_resi(1))%coo1(i,2)
			CA1(3)=Tgroup(ran_resi(1))%coo1(i,3)
		endif
	enddo

	call ran_gen(ran2,zero)
	ic=anint(ran2*180)
	delta_psi1=-(ic-90)
	cos_angle=cosd(delta_psi1)
	sin_angle=sind(delta_psi1)

	rotaxis_x=CA1(1)-C1(1)
	rotaxis_y=CA1(2)-C1(2)
	rotaxis_z=CA1(3)-C1(3)
	rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
	rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
	rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

	call axisrotation(rotaxis, cos_angle, sin_angle, m)

	do j=1, Tgroup(ran_resi(1))%cnum2
		Tgroup(ran_resi(1))%coo2(j,1)=Tgroup(ran_resi(1))%coo2(j,1)-C1(1)
		Tgroup(ran_resi(1))%coo2(j,2)=Tgroup(ran_resi(1))%coo2(j,2)-C1(2)
		Tgroup(ran_resi(1))%coo2(j,3)=Tgroup(ran_resi(1))%coo2(j,3)-C1(3)
	enddo
	do j=1, Tgroup(ran_resi(1))%cnum1
		Tgroup(ran_resi(1))%coo1(j,1)=Tgroup(ran_resi(1))%coo1(j,1)-C1(1)
		Tgroup(ran_resi(1))%coo1(j,2)=Tgroup(ran_resi(1))%coo1(j,2)-C1(2)
		Tgroup(ran_resi(1))%coo1(j,3)=Tgroup(ran_resi(1))%coo1(j,3)-C1(3)
	enddo
	Tgroupdata_backup(ran_resi(1))%coo(1)=Tgroupdata_backup(ran_resi(1))%coo(1)-C1(1)
	Tgroupdata_backup(ran_resi(1))%coo(2)=Tgroupdata_backup(ran_resi(1))%coo(2)-C1(2)
	Tgroupdata_backup(ran_resi(1))%coo(3)=Tgroupdata_backup(ran_resi(1))%coo(3)-C1(3)

	temp2=matmul(Tgroup(ran_resi(1))%coo2, m)
	Tgroup(ran_resi(1))%coo2=temp2

	temp1=matmul(Tgroup(ran_resi(1))%coo1, m)
	Tgroup(ran_resi(1))%coo1=temp1

	T_backup_temp=matmul(Tgroupdata_backup(ran_resi(1))%coo, m)
	Tgroupdata_backup(ran_resi(1))%coo=T_backup_temp

	do j=1, Tgroup(ran_resi(1))%cnum2
		Tgroup(ran_resi(1))%coo2(j,1)=anint((Tgroup(ran_resi(1))%coo2(j,1)+C1(1))*1000)/1000
		Tgroup(ran_resi(1))%coo2(j,2)=anint((Tgroup(ran_resi(1))%coo2(j,2)+C1(2))*1000)/1000
		Tgroup(ran_resi(1))%coo2(j,3)=anint((Tgroup(ran_resi(1))%coo2(j,3)+C1(3))*1000)/1000
	enddo
	do j=1, Tgroup(ran_resi(1))%cnum1
		Tgroup(ran_resi(1))%coo1(j,1)=anint((Tgroup(ran_resi(1))%coo1(j,1)+C1(1))*1000)/1000
		Tgroup(ran_resi(1))%coo1(j,2)=anint((Tgroup(ran_resi(1))%coo1(j,2)+C1(2))*1000)/1000
		Tgroup(ran_resi(1))%coo1(j,3)=anint((Tgroup(ran_resi(1))%coo1(j,3)+C1(3))*1000)/1000
	enddo
	Tgroupdata_backup(ran_resi(1))%coo(1)=anint((Tgroupdata_backup(ran_resi(1))%coo(1)+C1(1))*1000)/1000
	Tgroupdata_backup(ran_resi(1))%coo(2)=anint((Tgroupdata_backup(ran_resi(1))%coo(2)+C1(2))*1000)/1000
	Tgroupdata_backup(ran_resi(1))%coo(3)=anint((Tgroupdata_backup(ran_resi(1))%coo(3)+C1(3))*1000)/1000

	return
	end subroutine backbone_rotation_Nterm

	subroutine backbone_rotation_Cterm(group, groupdata_backup, ran_resi, delta_phi1, delta_psi1,  &
		delta_phi2, delta_psi2, Tgroup, Tgroupdata_backup)
	implicit none
	integer*8				:: ran_resi(3), i, j, k, l, temp_num, ic
	double precision					:: ran2
	double precision					:: N3(3), CA3(3), C3(3), N(3), CA(3), C(3)
	double precision					:: rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3)
	double precision					:: delta_phi1, delta_psi1, delta_phi2, delta_psi2, delta_phi3, delta_psi3
	double precision					:: delta_phi(2), delta_psi(2), cos_angle, sin_angle
	double precision					:: temp1(20,3), temp2(60,3), temp3(20,3), coo(20,3)
	double precision					:: T_backup_temp(3)
	character*4				:: atype(20)

	type(groupdetails)		:: group(gnum), Tgroup(gnum)
	type(databackup)		:: groupdata_backup(gnum), Tgroupdata_backup(gnum)

	delta_phi(1)=-delta_phi1; delta_phi(2)=-delta_phi2
	delta_psi(1)=-delta_psi1; delta_psi(2)=-delta_psi2

	Tgroup=group
	Tgroupdata_backup=groupdata_backup

	do i=1, 2
		cos_angle=cosd(delta_phi(i))
		sin_angle=sind(delta_phi(i))

		do k=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(k)=="N") then
				N(1)=Tgroup(ran_resi(i))%coo1(k,1)
				N(2)=Tgroup(ran_resi(i))%coo1(k,2)
				N(3)=Tgroup(ran_resi(i))%coo1(k,3)
			elseif(Tgroup(ran_resi(i))%atype1(k)=="CA") then
				CA(1)=Tgroup(ran_resi(i))%coo1(k,1)
				CA(2)=Tgroup(ran_resi(i))%coo1(k,2)
				CA(3)=Tgroup(ran_resi(i))%coo1(k,3)
			endif
		enddo

		rotaxis_x=CA(1)-N(1)
		rotaxis_y=CA(2)-N(2)
		rotaxis_z=CA(3)-N(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		temp_num=0
		do j=1, Tgroup(ran_resi(i))%cnum1
			if(Tgroup(ran_resi(i))%atype1(j)=="HA".or.Tgroup(ran_resi(i))%atype1(j)=="HA2" &
			.or.Tgroup(ran_resi(i))%atype1(j)=="HA3") then
				temp_num=temp_num+1
				atype(temp_num)=Tgroup(ran_resi(i))%atype1(j)
				coo(temp_num,1)=Tgroup(ran_resi(i))%coo1(j,1)-N(1)
				coo(temp_num,2)=Tgroup(ran_resi(i))%coo1(j,2)-N(2)
				coo(temp_num,3)=Tgroup(ran_resi(i))%coo1(j,3)-N(3)
			endif
		enddo

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=Tgroup(ran_resi(i))%coo2(j,1)-N(1)
			Tgroup(ran_resi(i))%coo2(j,2)=Tgroup(ran_resi(i))%coo2(j,2)-N(2)
			Tgroup(ran_resi(i))%coo2(j,3)=Tgroup(ran_resi(i))%coo2(j,3)-N(3)
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=Tgroup(ran_resi(i))%coo3(j,1)-N(1)
			Tgroup(ran_resi(i))%coo3(j,2)=Tgroup(ran_resi(i))%coo3(j,2)-N(2)
			Tgroup(ran_resi(i))%coo3(j,3)=Tgroup(ran_resi(i))%coo3(j,3)-N(3)
		enddo

		temp1=matmul(coo, m)
		coo=temp1

		temp2=matmul(Tgroup(ran_resi(i))%coo2, m)
		Tgroup(ran_resi(i))%coo2=temp2

		temp3=matmul(Tgroup(ran_resi(i))%coo3, m)
		Tgroup(ran_resi(i))%coo3=temp3

		do l=1, temp_num
			do j=1, Tgroup(ran_resi(i))%cnum1
				if(atype(l)==Tgroup(ran_resi(i))%atype1(j)) then
					Tgroup(ran_resi(i))%coo1(j,1)=anint((coo(l,1)+N(1))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,2)=anint((coo(l,2)+N(2))*1000)/1000
					Tgroup(ran_resi(i))%coo1(j,3)=anint((coo(l,3)+N(3))*1000)/1000
				endif
			enddo
		enddo

		do j=1, Tgroup(ran_resi(i))%cnum2
			Tgroup(ran_resi(i))%coo2(j,1)=anint((Tgroup(ran_resi(i))%coo2(j,1)+N(1))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,2)=anint((Tgroup(ran_resi(i))%coo2(j,2)+N(2))*1000)/1000
			Tgroup(ran_resi(i))%coo2(j,3)=anint((Tgroup(ran_resi(i))%coo2(j,3)+N(3))*1000)/1000
		enddo
		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=anint((Tgroup(ran_resi(i))%coo3(j,1)+N(1))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,2)=anint((Tgroup(ran_resi(i))%coo3(j,2)+N(2))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,3)=anint((Tgroup(ran_resi(i))%coo3(j,3)+N(3))*1000)/1000
		enddo

		do j=ran_resi(i)+1, ran_resi(3)
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-N(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-N(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-N(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-N(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-N(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-N(3)
			enddo
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-N(1)
				Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-N(2)
				Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-N(3)
			enddo
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-N(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-N(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-N(3)

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			temp3=matmul(Tgroup(j)%coo3, m)
			Tgroup(j)%coo3=temp3

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+N(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+N(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+N(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+N(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+N(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+N(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+N(1))*1000)/1000
				Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+N(2))*1000)/1000
				Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+N(3))*1000)/1000
			enddo
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+N(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+N(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+N(3))*1000)/1000
		enddo

		cos_angle=cosd(delta_psi(i))
		sin_angle=sind(delta_psi(i))

		do k=1, Tgroup(ran_resi(i))%cnum3
			if(Tgroup(ran_resi(i))%atype3(k)=="C") then
				C(1)=Tgroup(ran_resi(i))%coo3(k,1)
				C(2)=Tgroup(ran_resi(i))%coo3(k,2)
				C(3)=Tgroup(ran_resi(i))%coo3(k,3)
			endif
		enddo

		rotaxis_x=C(1)-CA(1)
		rotaxis_y=C(2)-CA(2)
		rotaxis_z=C(3)-CA(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=Tgroup(ran_resi(i))%coo3(j,1)-CA(1)
			Tgroup(ran_resi(i))%coo3(j,2)=Tgroup(ran_resi(i))%coo3(j,2)-CA(2)
			Tgroup(ran_resi(i))%coo3(j,3)=Tgroup(ran_resi(i))%coo3(j,3)-CA(3)
		enddo

		temp3=matmul(Tgroup(ran_resi(i))%coo3, m)
		Tgroup(ran_resi(i))%coo3=temp3

		do j=1, Tgroup(ran_resi(i))%cnum3
			Tgroup(ran_resi(i))%coo3(j,1)=anint((Tgroup(ran_resi(i))%coo3(j,1)+CA(1))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,2)=anint((Tgroup(ran_resi(i))%coo3(j,2)+CA(2))*1000)/1000
			Tgroup(ran_resi(i))%coo3(j,3)=anint((Tgroup(ran_resi(i))%coo3(j,3)+CA(3))*1000)/1000
		enddo

		do j=ran_resi(i)+1, ran_resi(3)
			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=Tgroup(j)%coo1(l,1)-CA(1)
				Tgroup(j)%coo1(l,2)=Tgroup(j)%coo1(l,2)-CA(2)
				Tgroup(j)%coo1(l,3)=Tgroup(j)%coo1(l,3)-CA(3)
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=Tgroup(j)%coo2(l,1)-CA(1)
				Tgroup(j)%coo2(l,2)=Tgroup(j)%coo2(l,2)-CA(2)
				Tgroup(j)%coo2(l,3)=Tgroup(j)%coo2(l,3)-CA(3)
			enddo
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=Tgroup(j)%coo3(l,1)-CA(1)
				Tgroup(j)%coo3(l,2)=Tgroup(j)%coo3(l,2)-CA(2)
				Tgroup(j)%coo3(l,3)=Tgroup(j)%coo3(l,3)-CA(3)
			enddo
			Tgroupdata_backup(j)%coo(1)=Tgroupdata_backup(j)%coo(1)-CA(1)
			Tgroupdata_backup(j)%coo(2)=Tgroupdata_backup(j)%coo(2)-CA(2)
			Tgroupdata_backup(j)%coo(3)=Tgroupdata_backup(j)%coo(3)-CA(3)

			temp1=matmul(Tgroup(j)%coo1, m)
			Tgroup(j)%coo1=temp1

			temp2=matmul(Tgroup(j)%coo2, m)
			Tgroup(j)%coo2=temp2

			temp3=matmul(Tgroup(j)%coo3, m)
			Tgroup(j)%coo3=temp3

			T_backup_temp=matmul(Tgroupdata_backup(j)%coo, m)
			Tgroupdata_backup(j)%coo=T_backup_temp

			do l=1, Tgroup(j)%cnum1
				Tgroup(j)%coo1(l,1)=anint((Tgroup(j)%coo1(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo1(l,2)=anint((Tgroup(j)%coo1(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo1(l,3)=anint((Tgroup(j)%coo1(l,3)+CA(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum2
				Tgroup(j)%coo2(l,1)=anint((Tgroup(j)%coo2(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo2(l,2)=anint((Tgroup(j)%coo2(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo2(l,3)=anint((Tgroup(j)%coo2(l,3)+CA(3))*1000)/1000
			enddo
			do l=1, Tgroup(j)%cnum3
				Tgroup(j)%coo3(l,1)=anint((Tgroup(j)%coo3(l,1)+CA(1))*1000)/1000
				Tgroup(j)%coo3(l,2)=anint((Tgroup(j)%coo3(l,2)+CA(2))*1000)/1000
				Tgroup(j)%coo3(l,3)=anint((Tgroup(j)%coo3(l,3)+CA(3))*1000)/1000
			enddo
			Tgroupdata_backup(j)%coo(1)=anint((Tgroupdata_backup(j)%coo(1)+CA(1))*1000)/1000
			Tgroupdata_backup(j)%coo(2)=anint((Tgroupdata_backup(j)%coo(2)+CA(2))*1000)/1000
			Tgroupdata_backup(j)%coo(3)=anint((Tgroupdata_backup(j)%coo(3)+CA(3))*1000)/1000
		enddo
	enddo

	do i=1, Tgroup(ran_resi(3))%cnum1
		if(Tgroup(ran_resi(3))%atype1(i)=="N") then
			N3(1)=Tgroup(ran_resi(3))%coo1(i,1)
			N3(2)=Tgroup(ran_resi(3))%coo1(i,2)
			N3(3)=Tgroup(ran_resi(3))%coo1(i,3)
		elseif(Tgroup(ran_resi(3))%atype1(i)=="CA") then
			CA3(1)=Tgroup(ran_resi(3))%coo1(i,1)
			CA3(2)=Tgroup(ran_resi(3))%coo1(i,2)
			CA3(3)=Tgroup(ran_resi(3))%coo1(i,3)
		endif
	enddo

	if(Tgroup(ran_resi(3))%gtype.ne."CPRO") then
		call ran_gen(ran2,zero)
		ic=anint(ran2*180)
		delta_phi3=-(ic-90)
		cos_angle=cosd(delta_phi3)
		sin_angle=sind(delta_phi3)

		rotaxis_x=CA3(1)-N3(1)
		rotaxis_y=CA3(2)-N3(2)
		rotaxis_z=CA3(3)-N3(3)
		rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

		call axisrotation(rotaxis, cos_angle, sin_angle, m)

		temp_num=0
		do j=1, Tgroup(ran_resi(3))%cnum1
			if(Tgroup(ran_resi(3))%atype1(j)=="HA".or.Tgroup(ran_resi(3))%atype1(j)=="HA2" & 
			.or.Tgroup(ran_resi(3))%atype1(j)=="HA3") then
				temp_num=temp_num+1
				atype(temp_num)=Tgroup(ran_resi(3))%atype1(j)
				coo(temp_num,1)=Tgroup(ran_resi(3))%coo1(j,1)-N3(1)
				coo(temp_num,2)=Tgroup(ran_resi(3))%coo1(j,2)-N3(2)
				coo(temp_num,3)=Tgroup(ran_resi(3))%coo1(j,3)-N3(3)
			endif
		enddo

		do j=1, Tgroup(ran_resi(3))%cnum2
			Tgroup(ran_resi(3))%coo2(j,1)=Tgroup(ran_resi(3))%coo2(j,1)-N3(1)
			Tgroup(ran_resi(3))%coo2(j,2)=Tgroup(ran_resi(3))%coo2(j,2)-N3(2)
			Tgroup(ran_resi(3))%coo2(j,3)=Tgroup(ran_resi(3))%coo2(j,3)-N3(3)
		enddo
		do j=1, Tgroup(ran_resi(3))%cnum3
			Tgroup(ran_resi(3))%coo3(j,1)=Tgroup(ran_resi(3))%coo3(j,1)-N3(1)
			Tgroup(ran_resi(3))%coo3(j,2)=Tgroup(ran_resi(3))%coo3(j,2)-N3(2)
			Tgroup(ran_resi(3))%coo3(j,3)=Tgroup(ran_resi(3))%coo3(j,3)-N3(3)
		enddo

		temp1=matmul(coo, m)
		coo=temp1

		temp2=matmul(Tgroup(ran_resi(3))%coo2, m)
		Tgroup(ran_resi(3))%coo2=temp2

		temp3=matmul(Tgroup(ran_resi(3))%coo3, m)
		Tgroup(ran_resi(3))%coo3=temp3

		do l=1, temp_num
			do j=1, Tgroup(ran_resi(3))%cnum1
				if(atype(l)==Tgroup(ran_resi(3))%atype1(j)) then
					Tgroup(ran_resi(3))%coo1(j,1)=anint((coo(l,1)+N3(1))*1000)/1000
					Tgroup(ran_resi(3))%coo1(j,2)=anint((coo(l,2)+N3(2))*1000)/1000
					Tgroup(ran_resi(3))%coo1(j,3)=anint((coo(l,3)+N3(3))*1000)/1000
				endif
			enddo
		enddo

		do j=1, Tgroup(ran_resi(3))%cnum2
			Tgroup(ran_resi(3))%coo2(j,1)=anint((Tgroup(ran_resi(3))%coo2(j,1)+N3(1))*1000)/1000
			Tgroup(ran_resi(3))%coo2(j,2)=anint((Tgroup(ran_resi(3))%coo2(j,2)+N3(2))*1000)/1000
			Tgroup(ran_resi(3))%coo2(j,3)=anint((Tgroup(ran_resi(3))%coo2(j,3)+N3(3))*1000)/1000
		enddo
		do j=1, Tgroup(ran_resi(3))%cnum3
			Tgroup(ran_resi(3))%coo3(j,1)=anint((Tgroup(ran_resi(3))%coo3(j,1)+N3(1))*1000)/1000
			Tgroup(ran_resi(3))%coo3(j,2)=anint((Tgroup(ran_resi(3))%coo3(j,2)+N3(2))*1000)/1000
			Tgroup(ran_resi(3))%coo3(j,3)=anint((Tgroup(ran_resi(3))%coo3(j,3)+N3(3))*1000)/1000
		enddo
	endif

	do i=1, Tgroup(ran_resi(3))%cnum3
		if(Tgroup(ran_resi(3))%atype3(i)=="C") then
			C3(1)=Tgroup(ran_resi(3))%coo3(i,1)
			C3(2)=Tgroup(ran_resi(3))%coo3(i,2)
			C3(3)=Tgroup(ran_resi(3))%coo3(i,3)
		endif
	enddo

	call ran_gen(ran2,zero)
	ic=anint(ran2*180)
	delta_psi3=-(ic-90)
	cos_angle=cosd(delta_psi3)
	sin_angle=sind(delta_psi3)

	rotaxis_x=C3(1)-CA3(1)
	rotaxis_y=C3(2)-CA3(2)
	rotaxis_z=C3(3)-CA3(3)
	rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
	rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
	rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)

	call axisrotation(rotaxis, cos_angle, sin_angle, m)

	do j=1, Tgroup(ran_resi(3))%cnum3
		Tgroup(ran_resi(3))%coo3(j,1)=Tgroup(ran_resi(3))%coo3(j,1)-CA3(1)
		Tgroup(ran_resi(3))%coo3(j,2)=Tgroup(ran_resi(3))%coo3(j,2)-CA3(2)
		Tgroup(ran_resi(3))%coo3(j,3)=Tgroup(ran_resi(3))%coo3(j,3)-CA3(3)
	enddo

	temp3=matmul(Tgroup(ran_resi(3))%coo3, m)
	Tgroup(ran_resi(3))%coo3=temp3

	do j=1, Tgroup(ran_resi(3))%cnum3
		Tgroup(ran_resi(3))%coo3(j,1)=anint((Tgroup(ran_resi(3))%coo3(j,1)+CA3(1))*1000)/1000
		Tgroup(ran_resi(3))%coo3(j,2)=anint((Tgroup(ran_resi(3))%coo3(j,2)+CA3(2))*1000)/1000
		Tgroup(ran_resi(3))%coo3(j,3)=anint((Tgroup(ran_resi(3))%coo3(j,3)+CA3(3))*1000)/1000
	enddo

	return
	end subroutine backbone_rotation_Cterm

	subroutine concertedrotation_center(group, groupdata_backup, ran_resi, phipsi_num, &
		group_candidates, groupdata_backup_candidates)
	implicit none
	integer*8						:: ran_resi(3), i, j, k, count, ic
	integer*8						:: phipsi_num, psi_num, H, feedback, flag
	double precision							:: C0(3), N1(3), CA1(3), C1(3), N2(3), CA2(3), C2(3) 
	double precision							:: N3(3), CA3(3), C3(3), N4(4), CA3_target(3)
	double precision							:: C10_temp(3), N1_temp(3), CA1_temp(3), C1_temp(3), N12_temp(3)
	double precision							:: C20_temp(3), N2_temp(3), CA2_temp(3), C2_temp(3), N22_temp(3)
	double precision							:: C30_temp(3), N3_temp(3), CA3_temp(3)
	double precision							:: c11, c22, c33, c44, dist, rmsd, ran2
	double precision							:: omg1, omg2, theta_omg0, theta_phi1, theta_psi1 
	double precision							:: theta_omg1, theta_phi2, theta_psi2, theta_omg2
	double precision							:: l_omg0, l_phi1, l_psi1, l_omg1, l_phi2, l_psi2, l_omg2, l_phi3
	double precision							:: U_phi(3), U_psi(3), U_omg(3), U_omg0(3), U_phi1(3), U_psi1(3)
	double precision							:: U_omg1(3), U_phi2(3), U_psi2(3), U_omg2(3), U_phi3(3)
	double precision							:: Lomg1(3), Lomg2(3), Q_omg1(3), Q_omg2(3), e(3), t(3) 
	double precision							:: t_prime(3), q2(3), q2_prime(3), m(3), n(3)
	double precision							:: Tlab(3,3), Tlab_inverse(3,3), Tphi1(3,3), Tpsi1(3,3), Tomg1(3,3)
	double precision							:: r_CA_N(3), r1(3), r1_old(3), r1_new(3), r_old_temp(3), r_new_temp(3)
	double precision							:: delta_theta, theta_old_temp, theta_new_temp
	double precision							:: delta_phi, phi1_old, phi1, phi2, phi(2), psi(2), psi1(2)
	double precision							:: phi1_temp, psi1_temp, phi2_temp, psi2_temp, phi3_temp, psi3_temp
	type(groupdetails)				:: group(gnum), group_candidates(20,gnum), Tgroup(gnum), Tgroup_new(gnum)
	type(databackup)				:: groupdata_backup(gnum), groupdata_backup_candidates(20, gnum) 
	type(databackup)                :: Tgroupdata_backup(gnum), Tgroupdata_backup_new(gnum)

	Tgroup=group
	Tgroupdata_backup=groupdata_backup
    
	phipsi_num=0

	do i=1, Tgroup(ran_resi(1))%cnum1
		if(Tgroup(ran_resi(1))%atype1(i)=="N") then
			N1(1)=Tgroup(ran_resi(1))%coo1(i,1)
			N1(2)=Tgroup(ran_resi(1))%coo1(i,2)
			N1(3)=Tgroup(ran_resi(1))%coo1(i,3)
		elseif(Tgroup(ran_resi(1))%atype1(i)=="CA") then
			CA1(1)=Tgroup(ran_resi(1))%coo1(i,1)
			CA1(2)=Tgroup(ran_resi(1))%coo1(i,2)
			CA1(3)=Tgroup(ran_resi(1))%coo1(i,3)
		endif
	enddo
	do i=1, Tgroup(ran_resi(1))%cnum3
		if(Tgroup(ran_resi(1))%atype3(i)=="C") then
			C1(1)=Tgroup(ran_resi(1))%coo3(i,1)
			C1(2)=Tgroup(ran_resi(1))%coo3(i,2)
			C1(3)=Tgroup(ran_resi(1))%coo3(i,3)
		endif
	enddo

	U_phi1(1)=CA1(1)-N1(1); U_phi1(2)=CA1(2)-N1(2); U_phi1(3)=CA1(3)-N1(3)
	U_psi1(1)=C1(1)-CA1(1); U_psi1(2)=C1(2)-CA1(2); U_psi1(3)=C1(3)-CA1(3)
	l_phi1=sqrt(U_phi1(1)**2+U_phi1(2)**2+U_phi1(3)**2)
	l_psi1=sqrt(U_psi1(1)**2+U_psi1(2)**2+U_psi1(3)**2)

	theta_phi1=acosd(dot_product(U_phi1, U_psi1)/(l_phi1*l_psi1))

	Tlab(1,1)=U_phi1(1)/l_phi1; Tlab(1,2)=U_phi1(2)/l_phi1; Tlab(1,3)=U_phi1(3)/l_phi1
	U_phi(1)=Tlab(1,1); U_phi(2)=Tlab(1,2); U_phi(3)=Tlab(1,3)

	U_psi(1)=U_psi1(1)/l_psi1; U_psi(2)=U_psi1(2)/l_psi1; U_psi(3)=U_psi1(3)/l_psi1

	Tlab(3,1)=(U_phi(2)*U_psi(3)-U_phi(3)*U_psi(2))/sind(theta_phi1)
	Tlab(3,2)=(U_phi(3)*U_psi(1)-U_phi(1)*U_psi(3))/sind(theta_phi1)
	Tlab(3,3)=(U_phi(1)*U_psi(2)-U_phi(2)*U_psi(1))/sind(theta_phi1)

	Tlab(2,1)=Tlab(3,2)*Tlab(1,3)-Tlab(3,3)*Tlab(1,2)
	Tlab(2,2)=Tlab(3,3)*Tlab(1,1)-Tlab(3,1)*Tlab(1,3)
	Tlab(2,3)=Tlab(3,1)*Tlab(1,2)-Tlab(3,2)*Tlab(1,1)

	Tlab_inverse=transpose(Tlab)

    !get angle between 60 and 74 (where delta_theta is between -4 and 4)
	count=0
	do while(count.le.20)
		call ran_gen(ran2,zero); delta_theta=ran2*8.0-4.0
		if((theta_phi1+delta_theta).ge.60.0.and.(theta_phi1+delta_theta).le.74.0) goto 30
		count=count+1
	enddo
	if(count.eq.21) goto 20
30	continue

	do ic=ran_resi(1), ran_resi(3)
		if(ic.ne.ran_resi(1)) then
			do i=1, Tgroup(ic)%cnum1
				r_old_temp(1)=Tgroup(ic)%coo1(i,1)-CA1(1)
				r_old_temp(2)=Tgroup(ic)%coo1(i,2)-CA1(2)
				r_old_temp(3)=Tgroup(ic)%coo1(i,3)-CA1(3)
				r1_old=matmul(Tlab, r_old_temp)
				theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
				if(r1_old(2).lt.0.0) then
					theta_old_temp=360.00-theta_old_temp
				endif
				theta_new_temp=theta_old_temp+delta_theta
				r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
				r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
				r1_new(3)=r1_old(3)
				r_new_temp=matmul(Tlab_inverse, r1_new)
				Tgroup(ic)%coo1(i,1)=r_new_temp(1)+CA1(1)
				Tgroup(ic)%coo1(i,2)=r_new_temp(2)+CA1(2)
				Tgroup(ic)%coo1(i,3)=r_new_temp(3)+CA1(3)
			enddo

			do i=1, Tgroup(ic)%cnum2
				r_old_temp(1)=Tgroup(ic)%coo2(i,1)-CA1(1)
				r_old_temp(2)=Tgroup(ic)%coo2(i,2)-CA1(2)
				r_old_temp(3)=Tgroup(ic)%coo2(i,3)-CA1(3)
				r1_old=matmul(Tlab, r_old_temp)
				theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
				if(r1_old(2).lt.0.0) then
					theta_old_temp=360.00-theta_old_temp
				endif
				theta_new_temp=theta_old_temp+delta_theta
				r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
				r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
				r1_new(3)=r1_old(3)
				r_new_temp=matmul(Tlab_inverse, r1_new)
				Tgroup(ic)%coo2(i,1)=r_new_temp(1)+CA1(1)
				Tgroup(ic)%coo2(i,2)=r_new_temp(2)+CA1(2)
				Tgroup(ic)%coo2(i,3)=r_new_temp(3)+CA1(3)
			enddo

			r_old_temp(1)=Tgroupdata_backup(ic)%coo(1)-CA1(1)
			r_old_temp(2)=Tgroupdata_backup(ic)%coo(2)-CA1(2)
			r_old_temp(3)=Tgroupdata_backup(ic)%coo(3)-CA1(3)
			r1_old=matmul(Tlab, r_old_temp)
			theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
			if(r1_old(2).lt.0.0) then
				theta_old_temp=360.00-theta_old_temp
			endif
			theta_new_temp=theta_old_temp+delta_theta
			r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
			r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
			r1_new(3)=r1_old(3)
			r_new_temp=matmul(Tlab_inverse, r1_new)
			Tgroupdata_backup(ic)%coo(1)=r_new_temp(1)+CA1(1)
			Tgroupdata_backup(ic)%coo(2)=r_new_temp(2)+CA1(2)
			Tgroupdata_backup(ic)%coo(3)=r_new_temp(3)+CA1(3)
		endif

		if(ic.ne.ran_resi(3)) then
			do i=1, Tgroup(ic)%cnum3
				r_old_temp(1)=Tgroup(ic)%coo3(i,1)-CA1(1)
				r_old_temp(2)=Tgroup(ic)%coo3(i,2)-CA1(2)
				r_old_temp(3)=Tgroup(ic)%coo3(i,3)-CA1(3)
				r1_old=matmul(Tlab, r_old_temp)
				theta_old_temp=acosd(r1_old(1)/sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2)))
				if(r1_old(2).lt.0.0) then
					theta_old_temp=360.00-theta_old_temp
				endif
				theta_new_temp=theta_old_temp+delta_theta
				r1_new(1)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*cosd(theta_new_temp)
				r1_new(2)=sqrt(r1_old(1)*r1_old(1)+r1_old(2)*r1_old(2))*sind(theta_new_temp)
				r1_new(3)=r1_old(3)
				r_new_temp=matmul(Tlab_inverse, r1_new)
				Tgroup(ic)%coo3(i,1)=r_new_temp(1)+CA1(1)
				Tgroup(ic)%coo3(i,2)=r_new_temp(2)+CA1(2)
				Tgroup(ic)%coo3(i,3)=r_new_temp(3)+CA1(3)
			enddo
		endif
	enddo

    !get position of C before N in first pivot
	do i=1, Tgroup(ran_resi(1)-1)%cnum3
		if(Tgroup(ran_resi(1)-1)%atype3(i)=="C") then
			C0(1)=Tgroup(ran_resi(1)-1)%coo3(i,1)
			C0(2)=Tgroup(ran_resi(1)-1)%coo3(i,2)
			C0(3)=Tgroup(ran_resi(1)-1)%coo3(i,3)
		endif
	enddo

    !get position of N and alpha carbon in first pivot
	do i=1, Tgroup(ran_resi(1))%cnum1
		if(Tgroup(ran_resi(1))%atype1(i)=="N") then
			N1(1)=Tgroup(ran_resi(1))%coo1(i,1)
			N1(2)=Tgroup(ran_resi(1))%coo1(i,2)
			N1(3)=Tgroup(ran_resi(1))%coo1(i,3)
		elseif(Tgroup(ran_resi(1))%atype1(i)=="CA") then
			CA1(1)=Tgroup(ran_resi(1))%coo1(i,1)
			CA1(2)=Tgroup(ran_resi(1))%coo1(i,2)
			CA1(3)=Tgroup(ran_resi(1))%coo1(i,3)
		endif
	enddo

    !get position of C in first pivot
	do i=1, Tgroup(ran_resi(1))%cnum3
		if(Tgroup(ran_resi(1))%atype3(i)=="C") then
			C1(1)=Tgroup(ran_resi(1))%coo3(i,1)
			C1(2)=Tgroup(ran_resi(1))%coo3(i,2)
			C1(3)=Tgroup(ran_resi(1))%coo3(i,3)
		endif
	enddo

    !get position of N and alpha carbon in second pivot
	do i=1, Tgroup(ran_resi(2))%cnum1
		if(Tgroup(ran_resi(2))%atype1(i)=="N") then
			N2(1)=Tgroup(ran_resi(2))%coo1(i,1)
			N2(2)=Tgroup(ran_resi(2))%coo1(i,2)
			N2(3)=Tgroup(ran_resi(2))%coo1(i,3)
		elseif(Tgroup(ran_resi(2))%atype1(i)=="CA") then
			CA2(1)=Tgroup(ran_resi(2))%coo1(i,1)
			CA2(2)=Tgroup(ran_resi(2))%coo1(i,2)
			CA2(3)=Tgroup(ran_resi(2))%coo1(i,3)
		endif
	enddo

    !get position of C in second pivot
	do i=1, Tgroup(ran_resi(2))%cnum3
		if(Tgroup(ran_resi(2))%atype3(i)=="C") then
			C2(1)=Tgroup(ran_resi(2))%coo3(i,1)
			C2(2)=Tgroup(ran_resi(2))%coo3(i,2)
			C2(3)=Tgroup(ran_resi(2))%coo3(i,3)
		endif
	enddo

    !get position of N and alpha carbon in third pivot
	do i=1, Tgroup(ran_resi(3))%cnum1
		if(Tgroup(ran_resi(3))%atype1(i)=="N") then
			N3(1)=Tgroup(ran_resi(3))%coo1(i,1)
			N3(2)=Tgroup(ran_resi(3))%coo1(i,2)
			N3(3)=Tgroup(ran_resi(3))%coo1(i,3)
		elseif(Tgroup(ran_resi(3))%atype1(i)=="CA") then
			CA3(1)=Tgroup(ran_resi(3))%coo1(i,1)
			CA3(2)=Tgroup(ran_resi(3))%coo1(i,2)
			CA3(3)=Tgroup(ran_resi(3))%coo1(i,3)
		endif
	enddo

    !get position of C in third pivot
	do i=1, Tgroup(ran_resi(3))%cnum3
		if(Tgroup(ran_resi(3))%atype3(i)=="C") then
			C3(1)=Tgroup(ran_resi(3))%coo3(i,1)
			C3(2)=Tgroup(ran_resi(3))%coo3(i,2)
			C3(3)=Tgroup(ran_resi(3))%coo3(i,3)
		endif
	enddo

    !get position of N in residue after third pivot
	do i=1, Tgroup(ran_resi(3)+1)%cnum1
		if(Tgroup(ran_resi(3)+1)%atype1(i)=="N") then
			N4(1)=Tgroup(ran_resi(3)+1)%coo1(i,1)
			N4(2)=Tgroup(ran_resi(3)+1)%coo1(i,2)
			N4(3)=Tgroup(ran_resi(3)+1)%coo1(i,3)
		endif
	enddo

    !calculate dihedral angles
	U_omg0(1)=N1(1)-C0(1); U_omg0(2)=N1(2)-C0(2); U_omg0(3)=N1(3)-C0(3)
	l_omg0=sqrt(U_omg0(1)**2+U_omg0(2)**2+U_omg0(3)**2)

	U_phi1(1)=CA1(1)-N1(1); U_phi1(2)=CA1(2)-N1(2); U_phi1(3)=CA1(3)-N1(3)
	U_psi1(1)=C1(1)-CA1(1); U_psi1(2)=C1(2)-CA1(2); U_psi1(3)=C1(3)-CA1(3)
	U_omg1(1)=N2(1)-C1(1); U_omg1(2)=N2(2)-C1(2); U_omg1(3)=N2(3)-C1(3)
	l_phi1=sqrt(U_phi1(1)**2+U_phi1(2)**2+U_phi1(3)**2)
	l_psi1=sqrt(U_psi1(1)**2+U_psi1(2)**2+U_psi1(3)**2)
	l_omg1=sqrt(U_omg1(1)**2+U_omg1(2)**2+U_omg1(3)**2)

	U_phi2(1)=CA2(1)-N2(1); U_phi2(2)=CA2(2)-N2(2); U_phi2(3)=CA2(3)-N2(3)
	U_psi2(1)=C2(1)-CA2(1); U_psi2(2)=C2(2)-CA2(2); U_psi2(3)=C2(3)-CA2(3)
	U_omg2(1)=N3(1)-C2(1); 	U_omg2(2)=N3(2)-C2(2); U_omg2(3)=N3(3)-C2(3)
	l_phi2=sqrt(U_phi2(1)**2+U_phi2(2)**2+U_phi2(3)**2)
	l_psi2=sqrt(U_psi2(1)**2+U_psi2(2)**2+U_psi2(3)**2)
	l_omg2=sqrt(U_omg2(1)**2+U_omg2(2)**2+U_omg2(3)**2)

	U_phi3(1)=CA3(1)-N3(1); U_phi3(2)=CA3(2)-N3(2); U_phi3(3)=CA3(3)-N3(3)
	l_phi3=sqrt(U_phi3(1)**2+U_phi3(2)**2+U_phi3(3)**2)

	call phipsiomg_angle(C0, N1, CA1, C1, phi1_old)
	call phipsiomg_angle(CA1, C1, N2, CA2, omg1)
	call phipsiomg_angle(CA2, C2, N3, CA3, omg2)

	theta_omg0=acosd(dot_product(U_omg0, U_phi1)/(l_omg0*l_phi1))
	theta_phi1=acosd(dot_product(U_phi1, U_psi1)/(l_phi1*l_psi1))
	theta_psi1=acosd(dot_product(U_psi1, U_omg1)/(l_psi1*l_omg1))
	theta_omg1=acosd(dot_product(U_omg1, U_phi2)/(l_omg1*l_phi2))
	theta_phi2=acosd(dot_product(U_phi2, U_psi2)/(l_phi2*l_psi2))
	theta_psi2=acosd(dot_product(U_psi2, U_omg2)/(l_psi2*l_omg2))
	theta_omg2=acosd(dot_product(U_omg2, U_phi3)/(l_omg2*l_phi3))

	Tlab(1,1)=U_phi1(1)/l_phi1; Tlab(1,2)=U_phi1(2)/l_phi1; Tlab(1,3)=U_phi1(3)/l_phi1
	U_phi(1)=Tlab(1,1); U_phi(2)=Tlab(1,2); U_phi(3)=Tlab(1,3)

	U_omg(1)=U_omg0(1)/l_omg0; U_omg(2)=U_omg0(2)/l_omg0; U_omg(3)=U_omg0(3)/l_omg0

	Tlab(3,1)=(U_phi(2)*U_omg(3)-U_phi(3)*U_omg(2))/sind(theta_omg0)
	Tlab(3,2)=(U_phi(3)*U_omg(1)-U_phi(1)*U_omg(3))/sind(theta_omg0)
	Tlab(3,3)=(U_phi(1)*U_omg(2)-U_phi(2)*U_omg(1))/sind(theta_omg0)

	Tlab(2,1)=Tlab(3,2)*Tlab(1,3)-Tlab(3,3)*Tlab(1,2)
	Tlab(2,2)=Tlab(3,3)*Tlab(1,1)-Tlab(3,1)*Tlab(1,3)
	Tlab(2,3)=Tlab(3,1)*Tlab(1,2)-Tlab(3,2)*Tlab(1,1)

	do i=1, group(ran_resi(3))%cnum1
		if(group(ran_resi(3))%atype1(i)=="CA") then
			CA3_target(1)=group(ran_resi(3))%coo1(i,1)
			CA3_target(2)=group(ran_resi(3))%coo1(i,2)
			CA3_target(3)=group(ran_resi(3))%coo1(i,3)
		endif
	enddo

	r_CA_N(1)=CA3_target(1)-N1(1); r_CA_N(2)=CA3_target(2)-N1(2); r_CA_N(3)=CA3_target(3)-N1(3)
	r1=matmul(Tlab,r_CA_N)
	r1(1)=r1(1)-l_phi1; r1(2)=r1(2); r1(3)=r1(3)

	Lomg1(1)=l_omg1+l_phi2*cosd(theta_omg1)
	Lomg1(2)=-l_phi2*sind(theta_omg1)*cosd(omg1)
	Lomg1(3)=l_phi2*sind(theta_omg1)*sind(omg1)

	Lomg2(1)=l_omg2+l_phi3*cosd(theta_omg2)
	Lomg2(2)=-l_phi3*sind(theta_omg2)*cosd(omg2)
	Lomg2(3)=l_phi3*sind(theta_omg2)*sind(omg2)

	Q_omg1(1)=l_psi1*cosd(theta_psi1)+Lomg1(1)
	Q_omg1(2)=l_psi1*sind(theta_psi1)+Lomg1(2)
	Q_omg1(3)=Lomg1(3)

	Q_omg2(1)=l_psi2*cosd(theta_psi2)+Lomg2(1)
	Q_omg2(2)=l_psi2*sind(theta_psi2)+Lomg2(2)
	Q_omg2(3)=Lomg2(3)

	e(1)=1; e(2)=0; e(3)=0

    !test different phi values to see if any lead to an acceptable backbone
	do delta_phi=delta_phi_min, delta_phi_max, delta_phi_step
		if(abs(delta_phi).le.2.0) goto 10
		phi1=DBLE(phi1_old+delta_phi)
		call transformatrix(theta_phi1, phi1, Tphi1)
		t=matmul(transpose(Tphi1), r1)

		t_prime(1)=(t(1)**2+t(2)**2+t(3)**2)+(Q_omg1(1)**2+Q_omg1(2)**2+Q_omg1(3)**2)-(Q_omg2(1)**2+Q_omg2(2)**2+ & 
					Q_omg2(3)**2)-2*t(1)*(cosd(theta_psi1)*Q_omg1(1)+sind(theta_psi1)*Q_omg1(2))

		t_prime(2)=2*t(2)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))+2*t(3)*Q_omg1(3)
		t_prime(3)=2*t(3)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))-2*t(2)*Q_omg1(3)

		c11=t_prime(1)/sqrt(t_prime(2)**2+t_prime(3)**2)
		if(abs(c11).le.1) then
			if(t_prime(3).gt.0) then
				H=0
			else
				H=180
			endif
			psi(1)=atand(t_prime(2)/t_prime(3))-asind(c11)+H
			psi(2)=atand(t_prime(2)/t_prime(3))+asind(c11)-180.0+H
			do j=1, 2
				do while (psi(j).le.(-179.5).or.psi(j).gt.(180.5))
					if(psi(j).le.(-179.5)) then
						psi(j)=psi(j)+360.0
					elseif(psi(j).gt.(180.5)) then
						psi(j)=psi(j)-360.0
					endif
				end do
			enddo

			psi_num=0
			do j=1, 2
				psi_num=psi_num+1
				psi1(psi_num)=psi(j)
			enddo
		else
			goto 10
		endif

		do i=1, psi_num
			call transformatrix(theta_psi1, psi1(i), Tpsi1)
			q2=matmul(transpose(Tpsi1), t)
			q2=q2-Q_omg1
			call transformatrix(theta_omg1, omg1, Tomg1)
			q2_prime=matmul(transpose(Tomg1), q2)

			m(1)=sind(theta_phi2)*(sind(theta_psi2)*Q_omg2(1)-cosd(theta_psi2)*Q_omg2(2))
			m(2)=sind(theta_phi2)*Q_omg2(3)
			m(3)=cosd(theta_phi2)*(cosd(theta_psi2)*Q_omg2(1)+sind(theta_psi2)*Q_omg2(2))-q2_prime(1)
			c22=m(3)/sqrt(m(1)**2+m(2)**2)
			if(abs(c22).le.1) then
				if(m(2).gt.0) then
					H=0
				else
					H=180
				endif
				psi(1)=atand(m(1)/m(2))-asind(c22)+H
				psi(2)=atand(m(1)/m(2))+asind(c22)-180.0+H
				do j=1, 2
					do while (psi(j).le.(-179.5).or.psi(j).gt.(180.5))
						if(psi(j).le.(-179.5)) then
							psi(j)=psi(j)+360.0
						elseif(psi(j).gt.(180.5)) then
							psi(j)=psi(j)-360.0
						endif
					end do
				enddo

				do j=1, 2
					n(1)=(sind(theta_phi2)*cosd(theta_psi2)+cosd(theta_phi2)*sind(theta_psi2)*cosd(psi(j)))*Q_omg2(1)
					n(1)=n(1)+(sind(theta_phi2)*sind(theta_psi2)-cosd(theta_phi2)*cosd(theta_psi2)*cosd(psi(j)))*Q_omg2(2)
					n(1)=n(1)-cosd(theta_phi2)*sind(psi(j))*Q_omg2(3)
					n(2)=sind(theta_psi2)*sind(psi(j))*Q_omg2(1)-cosd(theta_psi2)*sind(psi(j))*Q_omg2(2)+cosd(psi(j))*Q_omg2(3)
					n(3)=-q2_prime(2)
					c33=n(3)/sqrt(n(1)**2+n(2)**2)
					c44=q2_prime(3)
					if(abs(c33).le.1) then
						if(n(2).gt.0) then
							H=0
						else
							H=180
						endif
						phi(1)=atand(n(1)/n(2))-asind(c33)+H
						phi(2)=atand(n(1)/n(2))+asind(c33)-180.0+H
						do k=1, 2
							do while (phi(k).le.(-179.5).or.phi(k).gt.(180.5))
								if(phi(k).lt.(-179.5)) then
									phi(k)=phi(k)+360.0
								elseif(phi(k).gt.(180.5)) then
									phi(k)=phi(k)-360.0
								endif
							end do
						enddo

						flag=0
						do k=1, 2
							if(abs((n(1)*sind(phi(k))+n(2)*cosd(phi(k)))-c44).le.0.005) then
								phi2=phi(k)
								flag=1
							endif
						enddo

						if(flag==0) goto 40
						call backbone_rotation_center(Tgroup, Tgroupdata_backup, ran_resi, phi1, psi1(i), &
													phi2, psi(j), Tgroup_new, Tgroupdata_backup_new)

						do k=1, Tgroup_new(ran_resi(1)-1)%cnum3
							if(Tgroup_new(ran_resi(1)-1)%atype3(k)=="C") then
								C10_temp(1)=Tgroup_new(ran_resi(1)-1)%coo3(k,1)
								C10_temp(2)=Tgroup_new(ran_resi(1)-1)%coo3(k,2)
								C10_temp(3)=Tgroup_new(ran_resi(1)-1)%coo3(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(1))%cnum1
							if(Tgroup_new(ran_resi(1))%atype1(k)=="N") then
								N1_temp(1)=Tgroup_new(ran_resi(1))%coo1(k,1)
								N1_temp(2)=Tgroup_new(ran_resi(1))%coo1(k,2)
								N1_temp(3)=Tgroup_new(ran_resi(1))%coo1(k,3)
							elseif(Tgroup_new(ran_resi(1))%atype1(k)=="CA") then
								CA1_temp(1)=Tgroup_new(ran_resi(1))%coo1(k,1)
								CA1_temp(2)=Tgroup_new(ran_resi(1))%coo1(k,2)
								CA1_temp(3)=Tgroup_new(ran_resi(1))%coo1(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(1))%cnum3
							if(Tgroup_new(ran_resi(1))%atype3(k)=="C") then
								C1_temp(1)=Tgroup_new(ran_resi(1))%coo3(k,1)
								C1_temp(2)=Tgroup_new(ran_resi(1))%coo3(k,2)
								C1_temp(3)=Tgroup_new(ran_resi(1))%coo3(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(1)+1)%cnum1
							if(Tgroup_new(ran_resi(1)+1)%atype1(k)=="N") then
								N12_temp(1)=Tgroup_new(ran_resi(1)+1)%coo1(k,1)
								N12_temp(2)=Tgroup_new(ran_resi(1)+1)%coo1(k,2)
								N12_temp(3)=Tgroup_new(ran_resi(1)+1)%coo1(k,3)
							endif
						enddo

						do k=1, Tgroup_new(ran_resi(2)-1)%cnum3
							if(Tgroup_new(ran_resi(2)-1)%atype3(k)=="C") then
								C20_temp(1)=Tgroup_new(ran_resi(2)-1)%coo3(k,1)
								C20_temp(2)=Tgroup_new(ran_resi(2)-1)%coo3(k,2)
								C20_temp(3)=Tgroup_new(ran_resi(2)-1)%coo3(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(2))%cnum1
							if(Tgroup_new(ran_resi(2))%atype1(k)=="N") then
								N2_temp(1)=Tgroup_new(ran_resi(2))%coo1(k,1)
								N2_temp(2)=Tgroup_new(ran_resi(2))%coo1(k,2)
								N2_temp(3)=Tgroup_new(ran_resi(2))%coo1(k,3)
							elseif(Tgroup_new(ran_resi(2))%atype1(k)=="CA") then
								CA2_temp(1)=Tgroup_new(ran_resi(2))%coo1(k,1)
								CA2_temp(2)=Tgroup_new(ran_resi(2))%coo1(k,2)
								CA2_temp(3)=Tgroup_new(ran_resi(2))%coo1(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(2))%cnum3
							if(Tgroup_new(ran_resi(2))%atype3(k)=="C") then
								C2_temp(1)=Tgroup_new(ran_resi(2))%coo3(k,1)
								C2_temp(2)=Tgroup_new(ran_resi(2))%coo3(k,2)
								C2_temp(3)=Tgroup_new(ran_resi(2))%coo3(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(2)+1)%cnum1
							if(Tgroup_new(ran_resi(2)+1)%atype1(k)=="N") then
								N22_temp(1)=Tgroup_new(ran_resi(2)+1)%coo1(k,1)
								N22_temp(2)=Tgroup_new(ran_resi(2)+1)%coo1(k,2)
								N22_temp(3)=Tgroup_new(ran_resi(2)+1)%coo1(k,3)
							endif
						enddo

						do k=1, Tgroup_new(ran_resi(3)-1)%cnum3
							if(Tgroup_new(ran_resi(3)-1)%atype3(k)=="C") then
								C30_temp(1)=Tgroup_new(ran_resi(3)-1)%coo3(k,1)
								C30_temp(2)=Tgroup_new(ran_resi(3)-1)%coo3(k,2)
								C30_temp(3)=Tgroup_new(ran_resi(3)-1)%coo3(k,3)
							endif
						enddo
						do k=1, Tgroup_new(ran_resi(3))%cnum1
							if(Tgroup_new(ran_resi(3))%atype1(k)=="N") then
								N3_temp(1)=Tgroup_new(ran_resi(3))%coo1(k,1)
								N3_temp(2)=Tgroup_new(ran_resi(3))%coo1(k,2)
								N3_temp(3)=Tgroup_new(ran_resi(3))%coo1(k,3)
							elseif(Tgroup_new(ran_resi(3))%atype1(k)=="CA") then
								CA3_temp(1)=Tgroup_new(ran_resi(3))%coo1(k,1)
								CA3_temp(2)=Tgroup_new(ran_resi(3))%coo1(k,2)
								CA3_temp(3)=Tgroup_new(ran_resi(3))%coo1(k,3)
							endif
						enddo

                        !check if the new backbone bone angles satisfy the ramachandran map
						call phipsiomg_angle(C10_temp, N1_temp, CA1_temp, C1_temp, phi1_temp)
						call phipsiomg_angle(N1_temp, CA1_temp, C1_temp, N12_temp, psi1_temp)
						if(rama_map(anint(phi1_temp),anint(psi1_temp))==1.or.rama_map(anint(phi1_temp),anint(psi1_temp))==2) then
							call phipsiomg_angle(C20_temp, N2_temp, CA2_temp, C2_temp, phi2_temp)
							call phipsiomg_angle(N2_temp, CA2_temp, C2_temp, N22_temp, psi2_temp)
							if(rama_map(anint(phi2_temp),anint(psi2_temp))==1.or.rama_map(anint(phi2_temp),anint(psi2_temp))==2) then

                                !check the target positions are very close to the actual positions (i.e. residues that weren't supposed to move did not actually move)
								dist=sqrt((CA3_temp(1)-CA3_target(1))**2+(CA3_temp(2)-CA3_target(2))**2+(CA3_temp(3)-CA3_target(3))**2)
								if(dist.le.(0.05)) then
									call phipsiomg_angle(C30_temp, N3_temp, CA3_temp, C3, phi3_temp)
									call phipsiomg_angle(N3_temp, CA3_temp, C3, N4, psi3_temp)
									if(rama_map(anint(phi3_temp),anint(psi3_temp))==1.or.rama_map(anint(phi3_temp),anint(psi3_temp))==2) then
										call rmsd_calculation(Tgroup_new, rmsd)
										call backbonemove_criterion(rmsd, feedback) !check RMSD is within allowable range
										if(feedback==1) then
                                            !all checks passed, so add to possible list!
											phipsi_num=phipsi_num+1
											if(phipsi_num.gt.20) then
												phipsi_num=20
												goto 20
											endif
											do k=1, gnum
												group_candidates(phipsi_num,k)=Tgroup_new(k)
											enddo
											do k=1, gnum
												groupdata_backup_candidates(phipsi_num,k)%coo=Tgroupdata_backup_new(k)%coo
											enddo
										endif
									endif
								endif
							endif
						endif
					endif
40 					continue
				enddo
			endif
		enddo
10		continue

	enddo
20	continue

	return
	end subroutine concertedrotation_center

	subroutine concertedrotation_whole(group, groupdata_backup, ran_resi, phipsi_num, group_candidates, groupdata_backup_candidates)
	implicit none
	integer*8							:: ran_resi(3), k
	integer*8							:: phipsi_num, flag1, flag2

	type(groupdetails)				:: group(gnum), group_candidates(20,gnum)
	type(databackup)				:: groupdata_backup(gnum), groupdata_backup_candidates(20,gnum)

	type(groupdetails), dimension(:) :: Tgroup1(gnum), Tgroup2(gnum)
	type(databackup), dimension(:) :: Tgroupdata_backup1(gnum), Tgroupdata_backup2(gnum)

	phipsi_num=0
	call concertedrotation_N2Cterm(group, groupdata_backup, ran_resi, flag1, Tgroup1, Tgroupdata_backup1)
	if(flag1.ne.0) then
		call concertedrotation_C2Nterm(Tgroup1, Tgroupdata_backup1, ran_resi, flag2, Tgroup2, Tgroupdata_backup2)
		if(flag2.ne.0) then
			phipsi_num=phipsi_num+1
			do k=1, gnum
				group_candidates(phipsi_num,k)=Tgroup2(k)
			enddo
			do k=1, gnum
				groupdata_backup_candidates(phipsi_num,k)%coo=Tgroupdata_backup2(k)%coo
			enddo
		endif
	endif

	call concertedrotation_C2Nterm(group, groupdata_backup, ran_resi, flag1, Tgroup1, Tgroupdata_backup1)
	if(flag1.ne.0) then
		call concertedrotation_N2Cterm(Tgroup1, Tgroupdata_backup1, ran_resi, flag2, Tgroup2, Tgroupdata_backup2)
		if(flag2.ne.0) then
			phipsi_num=phipsi_num+1
			do k=1, gnum
				group_candidates(phipsi_num,k)=Tgroup2(k)
			enddo
			do k=1, gnum
				groupdata_backup_candidates(phipsi_num,k)%coo=Tgroupdata_backup2(k)%coo
			enddo
		endif
	endif

	return
	end subroutine concertedrotation_whole

	subroutine concertedrotation_Nterm(group, groupdata_backup, ran_resi, psiphi_num, group_candidates, groupdata_backup_candidates)
	implicit none
	integer*8							:: ran_resi(3), i, k
	integer*8							:: psiphi_num, feedback, Tresinum
	double precision							:: N30(3), C3(3), CA3(3), N3(3), C32(3)
	double precision							:: N20(3), C2(3), CA2(3), N2(3), C22(3)
	double precision							:: rmsd, ran2, ic, delta_psi3, delta_phi3, delta_psi2, delta_phi2
	double precision							:: psi3, phi3, psi2, phi2, Tpsi3, Tphi3, Tpsi2, Tphi2

	type(groupdetails)				:: group(gnum), group_candidates(20, gnum), Tgroup(gnum)
	type(databackup)				:: groupdata_backup(gnum), groupdata_backup_candidates(20,gnum), Tgroupdata_backup(gnum)

	do k=1, group(ran_resi(3)+1)%cnum1
		if(group(ran_resi(3)+1)%atype1(k)=="N") then
			N30(1)=group(ran_resi(3)+1)%coo1(k,1)
			N30(2)=group(ran_resi(3)+1)%coo1(k,2)
			N30(3)=group(ran_resi(3)+1)%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(3))%cnum3
		if(group(ran_resi(3))%atype3(k)=="C") then
			C3(1)=group(ran_resi(3))%coo3(k,1)
			C3(2)=group(ran_resi(3))%coo3(k,2)
			C3(3)=group(ran_resi(3))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(3))%cnum1
		if(group(ran_resi(3))%atype1(k)=="CA") then
			CA3(1)=group(ran_resi(3))%coo1(k,1)
			CA3(2)=group(ran_resi(3))%coo1(k,2)
			CA3(3)=group(ran_resi(3))%coo1(k,3)
		elseif(group(ran_resi(3))%atype1(k)=="N") then
			N3(1)=group(ran_resi(3))%coo1(k,1)
			N3(2)=group(ran_resi(3))%coo1(k,2)
			N3(3)=group(ran_resi(3))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(3)-1)%cnum3
		if(group(ran_resi(3)-1)%atype3(k)=="C") then
			C32(1)=group(ran_resi(3)-1)%coo3(k,1)
			C32(2)=group(ran_resi(3)-1)%coo3(k,2)
			C32(3)=group(ran_resi(3)-1)%coo3(k,3)
		endif
	enddo

	do k=1, group(ran_resi(2)+1)%cnum1
		if(group(ran_resi(2)+1)%atype1(k)=="N") then
			N20(1)=group(ran_resi(2)+1)%coo1(k,1)
			N20(2)=group(ran_resi(2)+1)%coo1(k,2)
			N20(3)=group(ran_resi(2)+1)%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum3
		if(group(ran_resi(2))%atype3(k)=="C") then
			C2(1)=group(ran_resi(2))%coo3(k,1)
			C2(2)=group(ran_resi(2))%coo3(k,2)
			C2(3)=group(ran_resi(2))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum1
		if(group(ran_resi(2))%atype1(k)=="CA") then
			CA2(1)=group(ran_resi(2))%coo1(k,1)
			CA2(2)=group(ran_resi(2))%coo1(k,2)
			CA2(3)=group(ran_resi(2))%coo1(k,3)
		elseif(group(ran_resi(2))%atype1(k)=="N") then
			N2(1)=group(ran_resi(2))%coo1(k,1)
			N2(2)=group(ran_resi(2))%coo1(k,2)
			N2(3)=group(ran_resi(2))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2)-1)%cnum3
		if(group(ran_resi(2)-1)%atype3(k)=="C") then
			C22(1)=group(ran_resi(2)-1)%coo3(k,1)
			C22(2)=group(ran_resi(2)-1)%coo3(k,2)
			C22(3)=group(ran_resi(2)-1)%coo3(k,3)
		endif
	enddo

	call phipsiomg_angle(N30, C3, CA3, N3, psi3)
	call phipsiomg_angle(C3, CA3, N3, C32, phi3)

	call phipsiomg_angle(N20, C2, CA2, N2, psi2)
	call phipsiomg_angle(C2, CA2, N2, C22, phi2)

	psiphi_num=0
	Tresinum=ran_resi(3)-ran_resi(1)
	do i=1, 20
		if(Tresinum.le.4) then
			call ran_gen(ran2,zero)
			ic=int(ran2*flavoredregion_number-1.0e-3)+1
			if(ic.gt.flavoredregion_number) ic=flavoredregion_number
			Tpsi3=flavored_region(ic,2)
			Tphi3=flavored_region(ic,1)

			delta_psi3=Tpsi3-psi3
			delta_phi3=Tphi3-phi3

			call ran_gen(ran2,zero)
			ic=int(ran2*flavoredregion_number-1.0e-3)+1
			if(ic.gt.flavoredregion_number) ic=flavoredregion_number
			Tpsi2=flavored_region(ic,2)
			Tphi2=flavored_region(ic,1)

			delta_psi2=Tpsi2-psi2
			delta_phi2=Tphi2-phi2
		else
			call ran_gen(ran2,zero); delta_psi3=ran2*6.9-3.45
			call ran_gen(ran2,zero); delta_phi3=ran2*6.9-3.45
			Tpsi3=psi3+delta_psi3; Tphi3=phi3+delta_phi3
			if(Tpsi3.ge.180.5) then
				Tpsi3=Tpsi3-360.0
				delta_psi3=Tpsi3-psi3
			elseif(Tpsi3.le.(-179.5)) then
				Tpsi3=Tpsi3+360.0
				delta_psi3=Tpsi3-psi3
			endif
			if(Tphi3.ge.180.5) then
				Tphi3=Tphi3-360.0
				delta_phi3=Tphi3-phi3
			elseif(Tphi3.le.(-179.5)) then
				Tphi3=Tphi3+360.0
				delta_phi3=Tphi3-phi3
			endif
			if(rama_map(anint(Tphi3), anint(Tpsi3))==0) goto 20

			call ran_gen(ran2,zero); delta_psi2=ran2*10.9-5.45
			call ran_gen(ran2,zero); delta_phi2=ran2*10.9-5.45
			Tpsi2=psi2+delta_psi2; Tphi2=phi2+delta_phi2
			if(Tpsi2.ge.180.5) then
				Tpsi2=Tpsi2-360.0
				delta_psi2=Tpsi2-psi2
			elseif(Tpsi2.le.(-179.5)) then
				Tpsi2=Tpsi2+360.0
				delta_psi2=Tpsi2-psi2
			endif
			if(Tphi2.ge.180.5) then
				Tphi2=Tphi2-360.0
				delta_phi2=Tphi2-phi2
			elseif(Tphi2.le.(-179.5)) then
				Tphi2=Tphi2+360.0
				delta_phi2=Tphi2-phi2
			endif
			if(rama_map(anint(Tphi2), anint(Tpsi2))==0) goto 20
		endif

		call backbone_rotation_Nterm(group, groupdata_backup, ran_resi, delta_psi3, delta_phi3, &
		delta_psi2, delta_phi2, Tgroup, Tgroupdata_backup)

		call rmsd_calculation(Tgroup, rmsd)
		call backbonemove_criterion(rmsd, feedback)
		if(feedback==1) then
			psiphi_num=psiphi_num+1
			do k=1, gnum
				group_candidates(psiphi_num,k)=Tgroup(k)
			enddo
			do k=1, gnum
				groupdata_backup_candidates(psiphi_num,k)%coo=Tgroupdata_backup(k)%coo
			enddo
			if(psiphi_num.eq.number4peptide_candidates) goto 10
		endif
20		continue
	enddo
10	continue

	return
	end subroutine concertedrotation_Nterm

	subroutine concertedrotation_C2Nterm(group, groupdata_backup, ran_resi, flag, Tgroup, Tgroupdata_backup)
	implicit none
	integer*8							:: ran_resi(3), i, k
	integer*8							:: flag, feedback
	double precision							:: N20(3), C2(3), CA2(3), N2(3), C22(3)
	double precision							:: rmsd, ran2, delta_psi3, delta_phi3, delta_psi2, delta_phi2
	double precision							:: psi2, phi2, Tpsi2, Tphi2

	type(groupdetails)				:: group(gnum), Tgroup(gnum)
	type(databackup)				:: groupdata_backup(gnum), Tgroupdata_backup(gnum)

	do k=1, group(ran_resi(2)+1)%cnum1
		if(group(ran_resi(2)+1)%atype1(k)=="N") then
			N20(1)=group(ran_resi(2)+1)%coo1(k,1)
			N20(2)=group(ran_resi(2)+1)%coo1(k,2)
			N20(3)=group(ran_resi(2)+1)%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum3
		if(group(ran_resi(2))%atype3(k)=="C") then
			C2(1)=group(ran_resi(2))%coo3(k,1)
			C2(2)=group(ran_resi(2))%coo3(k,2)
			C2(3)=group(ran_resi(2))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum1
		if(group(ran_resi(2))%atype1(k)=="CA") then
			CA2(1)=group(ran_resi(2))%coo1(k,1)
			CA2(2)=group(ran_resi(2))%coo1(k,2)
			CA2(3)=group(ran_resi(2))%coo1(k,3)
		elseif(group(ran_resi(2))%atype1(k)=="N") then
			N2(1)=group(ran_resi(2))%coo1(k,1)
			N2(2)=group(ran_resi(2))%coo1(k,2)
			N2(3)=group(ran_resi(2))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2)-1)%cnum3
		if(group(ran_resi(2)-1)%atype3(k)=="C") then
			C22(1)=group(ran_resi(2)-1)%coo3(k,1)
			C22(2)=group(ran_resi(2)-1)%coo3(k,2)
			C22(3)=group(ran_resi(2)-1)%coo3(k,3)
		endif
	enddo

	call phipsiomg_angle(N20, C2, CA2, N2, psi2)
	call phipsiomg_angle(C2, CA2, N2, C22, phi2)

	flag=0
	do i=1, 20
		delta_psi3=0.0
		call ran_gen(ran2,zero); delta_phi3=ran2*4.9-2.45

		call ran_gen(ran2,zero); delta_psi2=ran2*6.9-3.45
		call ran_gen(ran2,zero); delta_phi2=ran2*6.9-3.45
		Tpsi2=psi2+delta_psi2; Tphi2=phi2+delta_phi2
		if(Tpsi2.ge.180.5) then
			Tpsi2=Tpsi2-360.0
			delta_psi2=Tpsi2-psi2
		elseif(Tpsi2.le.(-179.5)) then
			Tpsi2=Tpsi2+360.0
			delta_psi2=Tpsi2-psi2
		endif
		if(Tphi2.ge.180.5) then
			Tphi2=Tphi2-360.0
			delta_phi2=Tphi2-phi2
		elseif(Tphi2.le.(-179.5)) then
			Tphi2=Tphi2+360.0
			delta_phi2=Tphi2-phi2
		endif
		if(rama_map(anint(Tphi2), anint(Tpsi2))==0) goto 20

		call backbone_rotation_Nterm(group, groupdata_backup, ran_resi, delta_psi3, delta_phi3, &
		delta_psi2, delta_phi2, Tgroup, Tgroupdata_backup)

		call rmsd_calculation(Tgroup, rmsd)
		call backbonemove_criterion(rmsd, feedback)
		if(feedback==1) then
			flag=flag+1
			goto 10
		endif
20		continue
	enddo
10	continue

	return
	end subroutine concertedrotation_C2Nterm

	subroutine concertedrotation_Cterm(group, groupdata_backup, ran_resi, phipsi_num, &
		group_candidates, groupdata_backup_candidates)
	implicit none
	integer*8							:: ran_resi(3), i, k
	integer*8							:: phipsi_num, feedback, Tresinum
	double precision							:: C10(3), N1(3), CA1(3), C1(3), N12(3)
	double precision							:: C20(3), N2(3), CA2(3), C2(3), N22(3)
	double precision							:: rmsd, ran2, ic, delta_phi1, delta_psi1, delta_phi2, delta_psi2
	double precision							:: phi1, psi1, phi2, psi2, Tphi1, Tpsi1, Tphi2, Tpsi2

	type(groupdetails)				:: group(gnum), group_candidates(20,gnum), Tgroup(gnum)
	type(databackup)				:: groupdata_backup(gnum), groupdata_backup_candidates(20,gnum), Tgroupdata_backup(gnum)

	do k=1, group(ran_resi(1)-1)%cnum3
		if(group(ran_resi(1)-1)%atype3(k)=="C") then
			C10(1)=group(ran_resi(1)-1)%coo3(k,1)
			C10(2)=group(ran_resi(1)-1)%coo3(k,2)
			C10(3)=group(ran_resi(1)-1)%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(1))%cnum1
		if(group(ran_resi(1))%atype1(k)=="N") then
			N1(1)=group(ran_resi(1))%coo1(k,1)
			N1(2)=group(ran_resi(1))%coo1(k,2)
			N1(3)=group(ran_resi(1))%coo1(k,3)
		elseif(group(ran_resi(1))%atype1(k)=="CA") then
			CA1(1)=group(ran_resi(1))%coo1(k,1)
			CA1(2)=group(ran_resi(1))%coo1(k,2)
			CA1(3)=group(ran_resi(1))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(1))%cnum3
		if(group(ran_resi(1))%atype3(k)=="C") then
			C1(1)=group(ran_resi(1))%coo3(k,1)
			C1(2)=group(ran_resi(1))%coo3(k,2)
			C1(3)=group(ran_resi(1))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(1)+1)%cnum1
		if(group(ran_resi(1)+1)%atype1(k)=="N") then
			N12(1)=group(ran_resi(1)+1)%coo1(k,1)
			N12(2)=group(ran_resi(1)+1)%coo1(k,2)
			N12(3)=group(ran_resi(1)+1)%coo1(k,3)
		endif
	enddo

	do k=1, group(ran_resi(2)-1)%cnum3
		if(group(ran_resi(2)-1)%atype3(k)=="C") then
			C20(1)=group(ran_resi(2)-1)%coo3(k,1)
			C20(2)=group(ran_resi(2)-1)%coo3(k,2)
			C20(3)=group(ran_resi(2)-1)%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum1
		if(group(ran_resi(2))%atype1(k)=="N") then
			N2(1)=group(ran_resi(2))%coo1(k,1)
			N2(2)=group(ran_resi(2))%coo1(k,2)
			N2(3)=group(ran_resi(2))%coo1(k,3)
		elseif(group(ran_resi(2))%atype1(k)=="CA") then
			CA2(1)=group(ran_resi(2))%coo1(k,1)
			CA2(2)=group(ran_resi(2))%coo1(k,2)
			CA2(3)=group(ran_resi(2))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum3
		if(group(ran_resi(2))%atype3(k)=="C") then
			C2(1)=group(ran_resi(2))%coo3(k,1)
			C2(2)=group(ran_resi(2))%coo3(k,2)
			C2(3)=group(ran_resi(2))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2)+1)%cnum1
		if(group(ran_resi(2)+1)%atype1(k)=="N") then
			N22(1)=group(ran_resi(2)+1)%coo1(k,1)
			N22(2)=group(ran_resi(2)+1)%coo1(k,2)
			N22(3)=group(ran_resi(2)+1)%coo1(k,3)
		endif
	enddo

	call phipsiomg_angle(C10, N1, CA1, C1, phi1)
	call phipsiomg_angle(N1, CA1, C1, N12, psi1)

	call phipsiomg_angle(C20, N2, CA2, C2, phi2)
	call phipsiomg_angle(N2, CA2, C2, N22, psi2)

	phipsi_num=0
	Tresinum=ran_resi(3)-ran_resi(1)
	do i=1, 20
		if(Tresinum.le.4) then
			call ran_gen(ran2,zero)
			ic=int(ran2*flavoredregion_number-1.0e-3)+1
			if(ic.gt.flavoredregion_number) ic=flavoredregion_number
			Tphi1=flavored_region(ic,1)
			Tpsi1=flavored_region(ic,2)

			delta_phi1=Tphi1-phi1
			delta_psi1=Tpsi1-psi1

			call ran_gen(ran2,zero)
			ic=int(ran2*flavoredregion_number-1.0e-3)+1
			if(ic.gt.flavoredregion_number) ic=flavoredregion_number
			Tphi2=flavored_region(ic,1)
			Tpsi2=flavored_region(ic,2)

			delta_phi2=Tphi2-phi2
			delta_psi2=Tpsi2-psi2
		else
			call ran_gen(ran2,zero); delta_phi1=ran2*6.8-3.4
			call ran_gen(ran2,zero); delta_psi1=ran2*6.8-3.4
			Tphi1=phi1+delta_phi1; Tpsi1=psi1+delta_psi1
			if(Tphi1.ge.180.5) then
				Tphi1=Tphi1-360.0
				delta_phi1=Tphi1-phi1
			elseif(Tphi1.le.(-179.5)) then
				Tphi1=Tphi1+360.0
				delta_phi1=Tphi1-phi1
			endif
			if(Tpsi1.ge.180.5) then
				Tpsi1=Tpsi1-360.0
				delta_psi1=Tpsi1-psi1
			elseif(Tpsi1.le.(-179.5)) then
				Tpsi1=Tpsi1+360.0
				delta_psi1=Tpsi1-psi1
			endif
			if(rama_map(anint(Tphi1), anint(Tpsi1))==0) goto 20

			call ran_gen(ran2,zero); delta_phi2=ran2*10.8-5.4
			call ran_gen(ran2,zero); delta_psi2=ran2*10.8-5.4
			Tphi2=phi2+delta_phi2; Tpsi2=psi2+delta_psi2
			if(Tphi2.ge.180.5) then
				Tphi2=Tphi2-360.0
				delta_phi2=Tphi2-phi2
			elseif(Tphi2.le.(-179.5)) then
				Tphi2=Tphi2+360.0
				delta_phi2=Tphi2-phi2
			endif
			if(Tpsi2.ge.180.5) then
				Tpsi2=Tpsi2-360.0
				delta_psi2=Tpsi2-psi2
			elseif(Tpsi2.le.(-179.5)) then
				Tpsi2=Tpsi2+360.0
				delta_psi2=Tpsi2-psi2
			endif
			if(rama_map(anint(Tphi2), anint(Tpsi2))==0) goto 20
		endif

		call backbone_rotation_Cterm(group, groupdata_backup, ran_resi, delta_phi1, delta_psi1,  & 
		delta_phi2, delta_psi2, Tgroup, Tgroupdata_backup)

		call rmsd_calculation(Tgroup, rmsd)
		call backbonemove_criterion(rmsd, feedback)
		if(feedback==1) then
			phipsi_num=phipsi_num+1
			do k=1, gnum
				group_candidates(phipsi_num,k)=Tgroup(k)
			enddo
			do k=1, gnum
				groupdata_backup_candidates(phipsi_num,k)%coo=Tgroupdata_backup(k)%coo
			enddo
			if(phipsi_num.eq.number4peptide_candidates) goto 10
		endif
20		continue
	enddo
10 	continue

	return
	end subroutine concertedrotation_Cterm

	subroutine concertedrotation_N2Cterm(group, groupdata_backup, ran_resi, flag, Tgroup, Tgroupdata_backup)
	implicit none
	integer*8							:: ran_resi(3), i, k
	integer*8							:: flag, feedback
	double precision							:: C20(3), N2(3), CA2(3), C2(3), N22(3)
	double precision							:: rmsd, ran2, delta_phi1, delta_psi1, delta_phi2, delta_psi2
	double precision							:: phi2, psi2, Tphi2, Tpsi2

	type(groupdetails)				:: group(gnum), Tgroup(gnum)
	type(databackup)				:: groupdata_backup(gnum), Tgroupdata_backup(gnum)

	do k=1, group(ran_resi(2)-1)%cnum3
		if(group(ran_resi(2)-1)%atype3(k)=="C") then
			C20(1)=group(ran_resi(2)-1)%coo3(k,1)
			C20(2)=group(ran_resi(2)-1)%coo3(k,2)
			C20(3)=group(ran_resi(2)-1)%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum1
		if(group(ran_resi(2))%atype1(k)=="N") then
			N2(1)=group(ran_resi(2))%coo1(k,1)
			N2(2)=group(ran_resi(2))%coo1(k,2)
			N2(3)=group(ran_resi(2))%coo1(k,3)
		elseif(group(ran_resi(2))%atype1(k)=="CA") then
			CA2(1)=group(ran_resi(2))%coo1(k,1)
			CA2(2)=group(ran_resi(2))%coo1(k,2)
			CA2(3)=group(ran_resi(2))%coo1(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2))%cnum3
		if(group(ran_resi(2))%atype3(k)=="C") then
			C2(1)=group(ran_resi(2))%coo3(k,1)
			C2(2)=group(ran_resi(2))%coo3(k,2)
			C2(3)=group(ran_resi(2))%coo3(k,3)
		endif
	enddo
	do k=1, group(ran_resi(2)+1)%cnum1
		if(group(ran_resi(2)+1)%atype1(k)=="N") then
			N22(1)=group(ran_resi(2)+1)%coo1(k,1)
			N22(2)=group(ran_resi(2)+1)%coo1(k,2)
			N22(3)=group(ran_resi(2)+1)%coo1(k,3)
		endif
	enddo

	call phipsiomg_angle(C20, N2, CA2, C2, phi2)
	call phipsiomg_angle(N2, CA2, C2, N22, psi2)

	flag=0
	do i=1, 20
		delta_phi1=0.0
		call ran_gen(ran2,zero); delta_psi1=ran2*4.9-2.45

		call ran_gen(ran2,zero); delta_phi2=ran2*6.9-3.45
		call ran_gen(ran2,zero); delta_psi2=ran2*6.9-3.45
		Tphi2=phi2+delta_phi2; Tpsi2=psi2+delta_psi2
		if(Tphi2.ge.180.5) then
			Tphi2=Tphi2-360.0
			delta_phi2=Tphi2-phi2
		elseif(Tphi2.le.(-179.5)) then
			Tphi2=Tphi2+360.0
			delta_phi2=Tphi2-phi2
		endif
		if(Tpsi2.ge.180.5) then
			Tpsi2=Tpsi2-360.0
			delta_psi2=Tpsi2-psi2
		elseif(Tpsi2.le.(-179.5)) then
			Tpsi2=Tpsi2+360.0
			delta_psi2=Tpsi2-psi2
		endif
		if(rama_map(anint(Tphi2), anint(Tpsi2))==0) goto 20

		call backbone_rotation_Cterm(group, groupdata_backup, ran_resi, delta_phi1, delta_psi1, &
		delta_phi2, delta_psi2, Tgroup, Tgroupdata_backup)

		call rmsd_calculation(Tgroup, rmsd)
		call backbonemove_criterion(rmsd, feedback)
		if(feedback==1) then
			flag=flag+1
			goto 10
		endif
20		continue
	enddo
10 	continue

	return
	end subroutine concertedrotation_N2Cterm

end module advancedfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module optimization_techniques

	use constant
	use randomgenerator
	use pdbfile
	use mathfunction
	use database
	use transplant
	use energy_calculation
	use advancedfunction
    use omp_lib
	use utilities

	contains
	subroutine MC_technique(energy_new, energy_old, vdw_new, vdw_old, ele_new, ele_old, sgb_new, sgb_old, &
	                        snp_new, snp_old, group, temp_group, feedback)
	implicit none
	integer*8							:: feedback
	double precision							:: ran2, energy_new, energy_old, energy_change
	double precision							:: vdw_old, ele_old, sgb_old, snp_old
	double precision							:: vdw_new, ele_new, sgb_new, snp_new
	type(groupdetails)				:: group(gnum), temp_group(gnum)

	feedback=0
	energy_change=energy_new-energy_old
	if(energy_change.le.0) then
		group=temp_group
		energy_old=energy_new
		vdw_old=vdw_new
		ele_old=ele_new
		sgb_old=sgb_new
		snp_old=snp_new
		feedback=1
	else
		call ran_gen(ran2,zero)
		if(ran2.le.exp(-energy_change/ekt_seq)) then
			group=temp_group
			energy_old=energy_new
			vdw_old=vdw_new
			ele_old=ele_new
			sgb_old=sgb_new
			snp_old=snp_new
			feedback=1
		endif
	endif

	return
	end subroutine MC_technique

	subroutine MC_technique_backbone(energy_new, energy_old, vdw_new, vdw_old, ele_new, ele_old, sgb_new, sgb_old, &
	                        snp_new, snp_old, group, temp_group, feedback)
	implicit none
	integer*8							:: feedback
	double precision							:: ran2, energy_new, energy_old, energy_change
	double precision							:: vdw_old, ele_old, sgb_old, snp_old
	double precision							:: vdw_new, ele_new, sgb_new, snp_new
	type(groupdetails)				:: group(:), temp_group(:)

	feedback=0
	energy_change=energy_new-energy_old
	if(energy_change.le.0) then
		group=temp_group
		energy_old=energy_new
		vdw_old=vdw_new
		ele_old=ele_new
		sgb_old=sgb_new
		snp_old=snp_new
		feedback=1
	else
		call ran_gen(ran2,zero)
		if(ran2.le.exp(-energy_change/ekt_backbone)) then
			group=temp_group
			energy_old=energy_new
			vdw_old=vdw_new
			ele_old=ele_new
			sgb_old=sgb_new
			snp_old=snp_new
			feedback=1
		endif
	endif

	return
	end subroutine MC_technique_backbone

	subroutine SCMF_technique(rotanum_1, rotanum_2, u, uij, i1_max, j2_max)
	implicit none
	integer*8			:: i, j, rotanum_1, rotanum_2
	integer*8			:: i1_max, j2_max
	double precision			:: u(2,40), uij(40,40), pmatrix_max
	double precision			:: pmatrix_old(2,40), pmatrix_new(2,40), effenergy(2,40), effenergy_total, error

	i1_max=0
	j2_max=0
	do i=1, rotanum_1
		pmatrix_old(1,i)=1.0/DBLE(rotanum_1)
	enddo
	do j=1, rotanum_2
		pmatrix_old(2,j)=1.0/DBLE(rotanum_2)
	enddo

10	continue
	pmatrix_new=pmatrix_old
	effenergy=0.0
	do i=1, rotanum_1
		effenergy(1,i)=u(1,i)
		do j=1, rotanum_2
			effenergy(1,i)=effenergy(1,i)+pmatrix_new(2,j)*uij(i,j)
		enddo
	enddo
	effenergy_total=0.0000001 !making slightly greater than 0 to prevent divide by 0 errors
	do i=1, rotanum_1
		effenergy_total=effenergy_total+exp(-effenergy(1,i)/ekt_scmf)
	enddo
	do i=1, rotanum_1
		pmatrix_new(1,i)=exp(-effenergy(1,i)/ekt_scmf)/effenergy_total
	enddo

	do j=1, rotanum_2
		effenergy(2,j)=u(2,j)
		do i=1, rotanum_1
			effenergy(2,j)=effenergy(2,j)+pmatrix_new(1,i)*uij(i,j)
		enddo
	enddo
	effenergy_total=0.0000001 !making slightly greater than 0 to prevent divide by 0 errors
	do j=1, rotanum_2
		effenergy_total=effenergy_total+exp(-effenergy(2,j)/ekt_scmf)
	enddo
	do j=1, rotanum_2
		pmatrix_new(2,j)=exp(-effenergy(2,j)/ekt_scmf)/effenergy_total
	enddo

	error=0.0
	do i=1, rotanum_1
		error=error+abs(pmatrix_new(1,i)-pmatrix_old(1,i))
	enddo
	do j=1, rotanum_2
		error=error+abs(pmatrix_new(2,j)-pmatrix_old(2,j))
	enddo

	if(error.gt.0.01) then
		pmatrix_old=lampda*pmatrix_new+(1-lampda)*pmatrix_old
		goto 10
	endif

	pmatrix_max=-1.0
	do i=1, rotanum_1
		if(pmatrix_max.le.pmatrix_new(1,i)) then
			pmatrix_max=pmatrix_new(1,i)
			i1_max=i
		endif
	enddo

	pmatrix_max=-1.0
	do j=1, rotanum_2
		if(pmatrix_max.le.pmatrix_new(2,j)) then
			pmatrix_max=pmatrix_new(2,j)
			j2_max=j
		endif
	enddo

	return
	end subroutine SCMF_technique

!below:
!switch can have several values, -1 tells it to always accept instead of choosing based on energy calculations. 0 tells it to accept/reject based on energy difference. from main function those are the only two options.


!from conrot section, switch can be 1-20, which refers to amino acid position taht was involved in the conrot change

	subroutine sequence_randomization(group, groupdata_backup)
	implicit none
	integer*8							:: attempt, i, j, rotanum_1, rotanum_2, feedback_1, feedback_2, flag1, flag2, flag3, trp_delta
	integer*8							:: ic1, ic2, ip, i1_max, j2_max, res_type_old, res_type_new
	integer*8							:: echou(2,40)
	double precision							:: ran2
	double precision							:: vdw_energy, vdw_energy_min
	double precision							:: u(2,40), uij(40,40)
	character*4						:: aminoacid_name_1, aminoacid_name_2
	character*4						:: group_name_1(3), group_name_2(3)

	type(groupdetails)				:: group(gnum)
	type(databackup)				:: groupdata_backup(gnum)
	type(groupdetails), dimension(:), allocatable :: temp_group_1, temp_group_2, aa_group_1, aa_group_2
	type(energyparameters), dimension(:), allocatable :: tgroup_para

	do attempt=1, sitenum  !!sitenum is equal to the number amino acids in the peptide that you wanna study
		call ran_gen(ran2,zero)
		if(ran2.le.scmfswitch) then 
			!we are going to replace one amino acid
5			continue
			call ran_gen(ran2,zero)
			ic1=int(ran2*sitenum-1.0e-3)+1 !choose an amino acid randomly
			if(ic1.gt.sitenum) ic1=sitenum
			
			call mc_choose_aminoacid(ic1, group, aminoacid_name_1, res_type_old, res_type_new, trp_delta) !this picks the residue which will replace the test residue
			allocate(aa_group_1(40))
			call findrotamer(ic1, group, aminoacid_name_1, rotanum_1, aa_group_1)  !find rotamers for the selected amino acid

			vdw_energy_min=1000.0	!initialize vdw to some high number

			ip=0 !this will store the index of the best rotamer
			allocate(temp_group_1(gnum))
			allocate(tgroup_para(gnum))
			do i=1, rotanum_1
				call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !do a trial of replacing the residue

				if(i==1) then
					call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters for the amino acid (only needs to happen once)
				endif

				call check_transplant(zero, ic1, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

				if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
	                call vdwenergy_sidechain_receptor(zero, ic1, zero, temp_group_1, tgroup_para, vdw_energy)

					if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
						vdw_energy_min=vdw_energy
						ip=i
					endif
				endif
			enddo

			if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule
				call residue_replace(ic1, group, groupdata_backup, ip, aa_group_1, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
				group=temp_group_1

				!since we changed the peptide, update the hydration properties and tryptophan count
				cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
				above_hyd(res_type_old) = above_hyd(res_type_old) - 1
				below_hyd(res_type_old) = below_hyd(res_type_old) + 1

				cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
				above_hyd(res_type_new) = above_hyd(res_type_new) + 1
				below_hyd(res_type_new) = below_hyd(res_type_new) - 1

				trp_count = trp_count + trp_delta
			endif
			deallocate(tgroup_para)
			deallocate(temp_group_1)
			deallocate(aa_group_1)

		else
			!we are going to swap two residues currently in the peptide
7			continue
			call ran_gen(ran2,zero)
			ic1=int(ran2*sitenum-1.0e-3)+1  !! pick 1 amino acid randomly
			if(ic1.gt.sitenum) ic1=sitenum
			ic2 = ic1
			do while(ic2.eq.ic1)
				!keep selecting residues until you pick one different from the first residue
				call ran_gen(ran2,zero)
				ic2=int(ran2*sitenum-1.0e-3)+1 
				if(ic2.gt.sitenum) ic2=sitenum
			enddo

			call groupinfo(group(ic1)%gtype, group_name_1, flag1)
			call groupinfo(group(ic2)%gtype, group_name_2, flag2)
			if(group_name_1(1)==group_name_2(1)) goto 7 !this tells the program to skip over the SCMF if the two residues have the same identity

			aminoacid_name_1=group_name_2(flag1) !load the name of the second residue into the slot for the first residue, and vice versa
			aminoacid_name_2=group_name_1(flag2)
			allocate(aa_group_1(40))
			allocate(aa_group_2(40))
			call findrotamer(ic1, group, aminoacid_name_1, rotanum_1, aa_group_1)!!calls rotamer possiblities
			call findrotamer(ic2, group, aminoacid_name_2, rotanum_2, aa_group_2)

			allocate(temp_group_1(gnum))
			allocate(temp_group_2(gnum))

			flag1=0 !
			allocate(tgroup_para(gnum))
			do i=1, rotanum_1
				call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !swap the residues

				if(i==1) then
					call energy_parameter(temp_group_1, tgroup_para, gnum) !load the force field parameters to the program, only do once since they are the same for all rotamers of a given amino acid
				endif

				call check_transplant(one, ic1, ic2, temp_group_1, feedback_1)!check for overlaps respond with feedback. if first argument =1, that means check residue site 1 with all other monomers except for ran_resi_2. has the same meaning in vdw energy below. this is because site 2 is currently a vacancy anyway

				if(feedback_1==1) then !if no atomic overlaps, then calculate vdw energy
					call vdwenergy_sidechain_receptor(one, ic1, ic2, temp_group_1, tgroup_para, vdw_energy)
					flag1=flag1+1
					u(1,flag1)=vdw_energy !assign the vdw energy value to u
					echou(1,flag1)=i !this tells index of the rotamer
				endif
			enddo
			deallocate(tgroup_para)

			if(flag1==0) goto 20

!! starting at flag2 below, this repeates the above steps for the 2nd residue
			flag2=0
			allocate(tgroup_para(gnum))
			do j=1, rotanum_2
				call residue_replace(ic2, group, groupdata_backup, j, aa_group_2, temp_group_2)

				if(j==1) then
					call energy_parameter(temp_group_2, tgroup_para, gnum)
				endif

				call check_transplant(one, ic2, ic1, temp_group_2, feedback_2)

				if(feedback_2==1) then
					call vdwenergy_sidechain_receptor(one, ic2, ic1, temp_group_2, tgroup_para, vdw_energy)
					flag2=flag2+1
					u(2,flag2)=vdw_energy
					echou(2,flag2)=j
				endif
			enddo
			deallocate(tgroup_para)

			if(flag2==0) goto 20
!!the above lines are for the 2nd residue, noted above in the comment about flag2

!flag3 loop is for caluclating energy among all possible rotamer combinations for both residues
			flag3=0
			allocate(tgroup_para(gnum))
			do i=1, flag1
				call residue_replace(ic1, group, groupdata_backup, echou(1,i), aa_group_1, temp_group_1)
				do j=1, flag2
					call residue_replace(ic2, temp_group_1, groupdata_backup, echou(2,j), aa_group_2, temp_group_2)

					if(i==1.and.j==1) then
						call energy_parameter(temp_group_2, tgroup_para, gnum)
					endif

					call check_transplant(zero, ic2, ic1, temp_group_2, feedback_2)

					if(feedback_2==1) then
						call vdwenergy_sidechain_receptor(one+1, ic2, ic1, temp_group_2, tgroup_para, vdw_energy)
						uij(i,j)=vdw_energy
						flag3=1
					else
						uij(i,j)=20.0
					endif
				enddo
			enddo
			deallocate(tgroup_para)

			if(flag3==0) goto 20

!this uses the SCMF to get the best combination, i.e. configuration with lowest vdw energy as weighted by rotamer probiabilities
			call SCMF_technique(flag1, flag2, u, uij, i1_max, j2_max)

			if(i1_max.ne.0.and.j2_max.ne.0) then
				call residue_replace(ic1, group, groupdata_backup, echou(1,i1_max), aa_group_1, temp_group_1)
				call residue_replace(ic2, temp_group_1, groupdata_backup, echou(2,j2_max), aa_group_2, temp_group_2)!!swap both residues according to the highest probability rotamer
			else
				goto 20
			endif

			call check_transplant(zero, ic1, ic2, temp_group_2, feedback_1) !check for atomic overlaps again for 1st residue
			if(feedback_1==0) then
				goto 20
			else
				call check_transplant(zero, ic2, ic1, temp_group_2, feedback_2) !check for atomic overlaps again for 2nd residue
				if(feedback_2==0) then
					goto 20
				endif
			endif

			group=temp_group_2

20			continue
			deallocate(aa_group_1)
			deallocate(aa_group_2)

			deallocate(temp_group_1)
			deallocate(temp_group_2)
		endif
	enddo
4	format(i7,i7,f18.13,a15)

	return
	end subroutine sequence_randomization

	subroutine sequence_optimization(group, SA_best_group, groupdata_backup, step, cycle, binding_energy_old, vdw_old, &
									ele_old, sgb_old, snp_old)
	implicit none
	integer*8						:: step, attempt, i, j, rotanum_1, rotanum_2, feedback_1, feedback_2, cycle, search_flag
	integer*8						:: ic1, ic2, ip, flag, flag1, flag2, stage, res_type_old, res_type_new, trp_delta
	integer*8						:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer*8						:: W_numex(atom_num), W_inb(atom_num,20), W_numex4(atom_num), W_inb4(atom_num, 60)
	double precision							:: ran2
	double precision							:: binding_energy_old, vdw_old, ele_old, sgb_old, snp_old
	double precision							:: binding_energy_new, vdw_new, ele_new, sgb_new, snp_new
	double precision							:: vdw_energy, vdw_energy_min, score, Tscore
	character*4						:: aminoacid_name_1, aminoacid_name_2
	character*4						:: group_name_1(3), group_name_2(3)
	double precision							:: ic1_prob, p_shift
	integer							:: pdb_flag, TABU_swap, TABU_pick1, TABU_pick2
	type(groupdetails)				:: group(gnum), Tgroup(gnum), SA_best_group(gnum)
	type(databackup)				:: groupdata_backup(gnum)
	type(groupdetails)				:: temp_group_1(gnum), temp_group_2(gnum), aa_group_1(40), aa_group_2(40)
	type(energyparameters)			:: tgroup_para(gnum)

	
	do attempt=1, sitenum  !!sitenum is equal to the number amino acids in the peptide
		call ran_gen(ran2,zero)
		if(ran2.le.scmfswitch) then !!scmfswitch is the value that chooses which way to go in the sequence changing part (change one amino acid?). this is able to be adjusted at the top of the code in the constant module.
5			continue
			!chose an amino acid to mutate based off res_selection vector 
			if (decomp_mutate_flag.eq.1) then
				call ran_gen(ran2,zero)
				search_flag = 0; ic1 = 1
				do while(search_flag.eq.0)
					if (ran2.gt.res_select_p_vector(ic1)) then
						ic1 = ic1+1	
					else
						search_flag = 1	
					endif
				enddo

			elseif (TABU_flag.eq.1) then
				call ran_gen(ran2,zero)
				TABU_pick1=min(int(ran2*TABU_allowed_length-1.0e-3)+1,sitenum) !pick residue that is currently allowed in TABU list
				ic1 = TABU_allowed_residues(TABU_pick1) !get the actual residue number this corresponds to			
			else
				call ran_gen(ran2,zero)
				ic1=min(int(ran2*sitenum-1.0e-3)+1,sitenum) !pick random residue
			endif

			SEQ_POINT_attempts = SEQ_POINT_attempts + 1
			SEQ_PerRes_attempts(ic1) = SEQ_PerRes_attempts(ic1) + 1

			call mc_choose_aminoacid(ic1, group, aminoacid_name_1, res_type_old, res_type_new, trp_delta) !this picks the residue which will replace the test residue
			call findrotamer(ic1, group, aminoacid_name_1, rotanum_1, aa_group_1)  !find how many rotamers for new amino acid

			vdw_energy_min=1000.0	!initialize
			ip=0 !initialize the rotamer index number
			if(aminoacid_name_1=="ALA".or.aminoacid_name_1=="PRO".or.aminoacid_name_1=="NALA".or.aminoacid_name_1=="NPRO".or. &
			aminoacid_name_1=="CALA".or.aminoacid_name_1=="CPRO".or.aminoacid_name_1=="GLY".or.aminoacid_name_1=="NGLY".or. &
			aminoacid_name_1=="CGLY") then
				do i=1, rotanum_1
					call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue

					if(i==1) then
						call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field parameters, the force field parameter doesnt change between rotamers so this only happens once
					endif

					call check_transplant(zero, ic1, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

					if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
						call vdwenergy_sidechain_receptor(zero, ic1, zero, temp_group_1, tgroup_para, vdw_energy)

						if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
							vdw_energy_min=vdw_energy
							ip=i
						endif
					endif
				enddo

				if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
					call residue_replace(ic1, group, groupdata_backup, ip, aa_group_1, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
					call atom_links(temp_group_1, W_numex, W_inb, W_numex4, W_inb4)
					call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
					                   binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
					call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
						              ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group, temp_group_1, feedback_2)

					if (feedback_2.eq.1) then
						!we actually changed the peptide, so update the hydration properties and tryptophan count
						cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
						above_hyd(res_type_old) = above_hyd(res_type_old) - 1
						below_hyd(res_type_old) = below_hyd(res_type_old) + 1

						cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
						above_hyd(res_type_new) = above_hyd(res_type_new) + 1
						below_hyd(res_type_new) = below_hyd(res_type_new) - 1

						trp_count = trp_count + trp_delta

						!check if top_score peptides array should be updated for ESURF calculation at end of program, and write pdb file if necessary
						if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
							call Update_Top_Peptides(binding_energy_old, cycle, step, attempt)
							call generatepdb(cycle, step, attempt, group)
						endif

						!update mutation attempt probabilities
						SEQ_POINT_success = SEQ_POINT_success + 1
						SEQ_PerRes_success(ic1) = SEQ_PerRes_success(ic1) + 1
						
						!update TABU search info
						if (TABU_flag.eq.1) then
							TABU_swap = TABU_allowed_residues(TABU_pick1)
							TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
							TABU_allowed_residues(TABU_allowed_length) = TABU_swap
							TABU_ban_lengths(TABU_pick1) = TABU_ban_lengths(TABU_allowed_length)
							TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
							TABU_allowed_length = TABU_allowed_length-1					
						endif

						if (binding_energy_old.lt.energy_min-.01) then
							energy_min = binding_energy_old
							if (SA_flag.eq.1) then
								!store group and hydration properties of current peptide
								SA_best_group = group
								SA_best_hyd = cur_hyd
								SA_above_hyd = above_hyd
								SA_below_hyd = below_hyd
								SA_trp_count = trp_count
							endif
							open(2, file="minimum_energy.txt", access="append")
								write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
							close(2)
							call generatepdb(cycle, step, attempt, SA_best_group)
						endif
					endif
				endif
			else
				flag1=0
				do i=1, rotanum_1
					call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue

					if(i==1) then
						call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
						call atom_links4sidechain(ic1, temp_group_1, S_numex, S_inb, S_numex4, S_inb4)
						call atom_links(temp_group_1, W_numex, W_inb, W_numex4, W_inb4)
					endif

					stage=0
					call sidechain_optimization(stage, ic1, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
					if(stage==1) then
!						call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
						if(flag1==0) then
							Tscore=score
							temp_group_2=temp_group_1
							flag1=1
						else
							if(score.le.Tscore) then
								Tscore=score
								temp_group_2=temp_group_1
							endif
						endif
					endif
				enddo

				if(flag1==1) then
					call sidechain_optimization4binding(ic1, temp_group_2, tgroup_para, W_numex, W_inb, W_numex4, W_inb4)
					call check_transplant(zero, ic1, zero, temp_group_2, feedback_1)
					if(feedback_1==1) then
						call bindingenergy(temp_group_2, tgroup_para, W_numex, W_inb, W_numex4, W_inb4,  &
						binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
						call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
										  ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group, temp_group_2, feedback_2)

						if (feedback_2.eq.1) then
							!we actually changed the peptide, so update the hydration properties
							cur_hyd(res_type_old) = cur_hyd(res_type_old) - 1
							above_hyd(res_type_old) = above_hyd(res_type_old) - 1
							below_hyd(res_type_old) = below_hyd(res_type_old) + 1
	
							cur_hyd(res_type_new) = cur_hyd(res_type_new) + 1
							above_hyd(res_type_new) = above_hyd(res_type_new) + 1
							below_hyd(res_type_new) = below_hyd(res_type_new) - 1

							trp_count = trp_count + trp_delta

							!check if top_score peptides array should be updated for ESURF calculation at end of program, and write pdb file if necessary
							if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
								call Update_Top_Peptides(binding_energy_old, cycle, step, attempt)
								call generatepdb(cycle, step, attempt, group)
							endif

							!update mutation attempt probabilities
							SEQ_POINT_success = SEQ_POINT_success + 1
							SEQ_PerRes_success(ic1) = SEQ_PerRes_success(ic1) + 1

							!update TABU search info
							if (TABU_flag.eq.1) then
								TABU_swap = TABU_allowed_residues(TABU_pick1)
								TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
								TABU_allowed_residues(TABU_allowed_length) = TABU_swap
								TABU_ban_lengths(TABU_pick1) = TABU_ban_lengths(TABU_allowed_length)
								TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
								TABU_allowed_length = TABU_allowed_length-1					
							endif

							!check if this the lowest energy we have found, and update simulated annealing best group (if applicable)
							if (binding_energy_old.lt.energy_min-.01) then
								energy_min = binding_energy_old
								if (SA_flag.eq.1) then
									SA_best_group = group
									SA_best_hyd = cur_hyd
									SA_above_hyd = above_hyd
									SA_below_hyd = below_hyd
									SA_trp_count = trp_count
								endif
								open(2, file="minimum_energy.txt", access="append")
									write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
								close(2)
								call generatepdb(cycle, step, attempt, SA_best_group)
							endif
						endif
					endif
				endif
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "Not SCMF" !!not SCMF indicates that it was a one amino change instead of swapping amino acids
				!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
				write(3,*) (group(j)%gtype, j=1, sitenum)
				write(3,*) "*******************************"
			close(3)

		else
			!we are going to swap two amino acids already in the peptide
			if (decomp_mutate_flag.eq.1) then
				!select the first residue per probability vector
				call ran_gen(ran2,zero)
				search_flag = 0; ic1 = 1
				do while(search_flag.eq.0)
					if (ran2.gt.res_select_p_vector(ic1)) then
						ic1 = ic1+1	
					else
						search_flag = 1	
					endif
				enddo

				!select second residue per probability vector (ensuring second residue is different from first residue)
				!do this in a way such that residue 1 is not reselected and are guaranteed to select a new residue 
				call ran_gen(ran2,zero)
				search_flag = 0; ic2 = 1; 
				if (ic1.eq.1) then
					ic1_prob = res_select_p_vector(ic1)
				else
					ic1_prob = res_select_p_vector(ic1) - res_select_p_vector(ic1-1)
				endif
				ran2 = ran2*(1 - ic1_prob) !this shifts cap on probability so we always are able to select a second residue
				p_shift = 0.0
				do while(search_flag.eq.0)
					if (ran2.gt.(res_select_p_vector(ic2)-p_shift).and.ic2.ne.ic1) then
						ic2 = ic2 + 1
						if (ic2.eq.ic1) then
							ic2 = ic2 + 1	
							p_shift = ic1_prob
						endif
					else
						search_flag = 1	
					endif
				enddo

			elseif (TABU_flag.eq.1) then
				call ran_gen(ran2,zero)
				TABU_pick1=min(int(ran2*TABU_allowed_length-1.0e-3)+1,TABU_allowed_length)   !! ic1 is randomly picked first
				do while(.true.)
					call ran_gen(ran2,zero)
					TABU_pick2=min(int(ran2*TABU_allowed_length-1.0e-3)+1,TABU_allowed_length)!!ic2 is the second amino acid selected for exchange, randomly chosen
					if(TABU_pick1.ne.TABU_pick2) then !!this makes sure that you don't pick the same residue twice
						exit
					endif
				enddo
				ic1 = TABU_allowed_residues(TABU_pick1)
				ic2 = TABU_allowed_residues(TABU_pick2)

			else 
				call ran_gen(ran2,zero)
				ic1=int(ran2*sitenum-1.0e-3)+1  !! ic1 is randomly picked first
				if(ic1.gt.sitenum) ic1=sitenum
				do while(.true.)
					call ran_gen(ran2,zero)
					ic2=int(ran2*sitenum-1.0e-3)+1 !!ic2 is the second amino acid selected for exchange, randomly chosen
					if(ic2.gt.sitenum) ic2=sitenum
					if(ic1.ne.ic2) then !!this makes sure that you don't pick the same residue twice
						exit
					endif
				enddo
			endif

			SEQ_SCMF_attempts = SEQ_SCMF_attempts + 1
			SEQ_PerRes_attempts(ic1) = SEQ_PerRes_attempts(ic1) + 1
			SEQ_PerRes_attempts(ic2) = SEQ_PerRes_attempts(ic2) + 1

			call groupinfo(group(ic1)%gtype, group_name_1, flag1)
			call groupinfo(group(ic2)%gtype, group_name_2, flag2)

			aminoacid_name_1=group_name_2(flag1) !load the name of the second residue into the slot for the first residue, and vice versa
			aminoacid_name_2=group_name_1(flag2)
			call findrotamer(ic1, group, aminoacid_name_1, rotanum_1, aa_group_1)!!calls rotamer possiblities
			call findrotamer(ic2, group, aminoacid_name_2, rotanum_2, aa_group_2)

			flag=0
			vdw_energy_min=1000.0	!initialize
			ip=0 !initialize the rotamer index number
			if(aminoacid_name_1=="GLY".or.aminoacid_name_1=="ALA".or.aminoacid_name_1=="PRO".or. &
			aminoacid_name_1=="NGLY".or.aminoacid_name_1=="NALA".or.aminoacid_name_1=="NPRO".or.&
			aminoacid_name_1=="CGLY".or.aminoacid_name_1=="CALA".or.aminoacid_name_1=="CPRO") then
				do i=1, rotanum_1
					call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue

					if(i==1) then
						call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
					endif

					call check_transplant(zero, ic1, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

					if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
						call vdwenergy_sidechain_receptor(zero, ic1, zero, temp_group_1, tgroup_para, vdw_energy)

						if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
							vdw_energy_min=vdw_energy
							ip=i
						endif
					endif
				enddo

				if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
					call residue_replace(ic1, group, groupdata_backup, ip, aa_group_1, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
					Tgroup=temp_group_1
					flag=1
				endif
			else
				flag1=0
				do i=1, rotanum_1
					call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue

					if(i==1) then
						call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
						call atom_links4sidechain(ic1, temp_group_1, S_numex, S_inb, S_numex4, S_inb4)
						call atom_links(temp_group_1, W_numex, W_inb, W_numex4, W_inb4)
					endif

					stage=0
					call sidechain_optimization(stage, ic1, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
					if(stage==1) then
!						call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
						if(flag1==0) then
							Tscore=score
							temp_group_2=temp_group_1
							flag1=1
						else
							if(score.le.Tscore) then
								Tscore=score
								temp_group_2=temp_group_1
							endif
						endif
					endif
				enddo

				if(flag1==1) then
					call sidechain_optimization4binding(ic1, temp_group_2, tgroup_para, W_numex, W_inb, W_numex4, W_inb4)
					call check_transplant(zero, ic1, zero, temp_group_2, feedback_1)
					if(feedback_1==1) then
						Tgroup=temp_group_2
						flag=1
					endif
				endif
			endif

			if(flag==1) then
				vdw_energy_min=1000.0	!initialize
				ip=0 !initialize the rotamer index number
				if(aminoacid_name_2=="GLY".or.aminoacid_name_2=="ALA".or.aminoacid_name_2=="PRO" &
				.or.aminoacid_name_2=="NGLY".or.aminoacid_name_2=="NALA".or.aminoacid_name_2=="NPRO".or. & 
				aminoacid_name_2=="CGLY".or.aminoacid_name_2=="CALA".or.aminoacid_name_2=="CPRO") then
					do i=1, rotanum_2
						call residue_replace(ic2, Tgroup, groupdata_backup, i, aa_group_2, temp_group_1) !replaces the residue

						if(i==1) then
							call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
						endif

						call check_transplant(zero, ic2, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

						if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
							call vdwenergy_sidechain_receptor(zero, ic2, zero, temp_group_1, tgroup_para, vdw_energy)

							if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
								vdw_energy_min=vdw_energy
								ip=i
							endif
						endif
					enddo

					if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
						call residue_replace(ic2, Tgroup, groupdata_backup, ip, aa_group_2, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
						call atom_links(temp_group_1, W_numex, W_inb, W_numex4, W_inb4)
						call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
										   binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
						call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
										  ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group, temp_group_1, feedback_2)

						!check if top_score peptides array should be updated for ESURF calculation at end of program, and write pdb file if necessary
						if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
							call Update_Top_Peptides(binding_energy_old, cycle, step, attempt)
							call generatepdb(cycle, step, attempt, group)
						endif

						!update mutation attempt probabilities
						if (feedback_2.eq.1) then
							SEQ_SCMF_success = SEQ_SCMF_success + 1
							SEQ_PerRes_success(ic1) = SEQ_PerRes_success(ic1) + 1
							SEQ_PerRes_success(ic2) = SEQ_PerRes_success(ic2) + 1
							
							!update TABU search info
							if (TABU_flag.eq.1) then
								if (TABU_flag.eq.1) then
									TABU_swap = TABU_allowed_residues(TABU_pick1)
									TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
									TABU_allowed_residues(TABU_allowed_length) = TABU_swap
									TABU_swap = TABU_allowed_residues(TABU_pick2)
									TABU_allowed_residues(TABU_pick2) = TABU_allowed_residues(TABU_allowed_length-1)
									TABU_allowed_residues(TABU_allowed_length-1) = TABU_swap
									TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
									TABU_ban_lengths(TABU_allowed_length-1) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
									TABU_allowed_length = TABU_allowed_length-2	
								endif				
							endif
						endif
						
						!check if best peptide score should be updated
						if (binding_energy_old.lt.energy_min-.01) then
							energy_min = binding_energy_old
							if (SA_flag.eq.1) then
								SA_best_group = group
								SA_best_hyd = cur_hyd
								SA_above_hyd = above_hyd
								SA_below_hyd = below_hyd
								SA_trp_count = trp_count
							endif
							open(2, file="minimum_energy.txt", access="append")
								write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
							close(2)
							call generatepdb(cycle, step, attempt, SA_best_group)
						endif
					endif
				else
					flag1=0
					do i=1, rotanum_2
						call residue_replace(ic2, Tgroup, groupdata_backup, i, aa_group_2, temp_group_1) !replaces the residue

						if(i==1) then
							call energy_parameter(temp_group_1, tgroup_para, gnum) !loads the force field paramters, the force field parameter doesnt change between rotamers so this only happens once
							call atom_links4sidechain(ic2, temp_group_1, S_numex, S_inb, S_numex4, S_inb4)
							call atom_links(temp_group_1, W_numex, W_inb, W_numex4, W_inb4)
						endif

						stage=0
						call sidechain_optimization(stage, ic2, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
						if(stage==1) then
!							call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
							if(flag1==0) then
								Tscore=score
								temp_group_2=temp_group_1
								flag1=1
							else
								if(score.le.Tscore) then
									Tscore=score
									temp_group_2=temp_group_1
								endif
							endif
						endif
					enddo

					if(flag1==1) then
						call sidechain_optimization4binding(ic2, temp_group_2, tgroup_para, W_numex, W_inb, W_numex4, W_inb4)
						call check_transplant(zero, ic2, zero, temp_group_2, feedback_1)
						if(feedback_1==1) then
							call bindingenergy(temp_group_2, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
							binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
							call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
											  ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group, temp_group_2, feedback_2)
							
							!check if top_score peptides array should be updated for ESURF calculation at end of program, and write pdb file if necessary
							if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
								call Update_Top_Peptides(binding_energy_old, cycle, step, attempt)
								call generatepdb(cycle, step, attempt, group)
							endif

							!update mutation attempt probabilities
							if (feedback_2.eq.1) then
								SEQ_SCMF_success = SEQ_SCMF_success + 1
								SEQ_PerRes_success(ic1) = SEQ_PerRes_success(ic1) + 1
								SEQ_PerRes_success(ic2) = SEQ_PerRes_success(ic2) + 1
								
								!update TABU search info
								if (TABU_flag.eq.1) then
									TABU_swap = TABU_allowed_residues(TABU_pick1)
									TABU_allowed_residues(TABU_pick1) = TABU_allowed_residues(TABU_allowed_length)
									TABU_allowed_residues(TABU_allowed_length) = TABU_swap
									TABU_swap = TABU_allowed_residues(TABU_pick2)
									TABU_allowed_residues(TABU_pick2) = TABU_allowed_residues(TABU_allowed_length-1)
									TABU_allowed_residues(TABU_allowed_length-1) = TABU_swap
									TABU_ban_lengths(TABU_allowed_length) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
									TABU_ban_lengths(TABU_allowed_length-1) = TABU_ban_length+1 !add 1 extra, since all values in array decreased by 1 before going back through loop
									TABU_allowed_length = TABU_allowed_length-2	
								endif
							endif

							if (binding_energy_old.lt.energy_min-.01) then
								energy_min = binding_energy_old
								if (SA_flag.eq.1) then
									SA_best_group = group
									SA_best_hyd = cur_hyd
									SA_above_hyd = above_hyd
									SA_below_hyd = below_hyd
									SA_trp_count = trp_count
								endif
								open(2, file="minimum_energy.txt", access="append")
									write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
								close(2)
								call generatepdb(cycle, step, attempt, SA_best_group)
							endif
						endif
					endif
				endif
			endif

			open(3, file="energydetails.txt",  access="append")!!output the energy information
				write(3,4) cycle,  step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "SCMF"
				!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
				write(3,*) (group(j)%gtype, j=1, sitenum)
				write(3,*) "*******************************"
			close(3)
		endif
		if (TABU_flag.eq.1.and.feedback_2.eq.1) call update_TABU_ban_list() !only update TABU list whenever a sequence change is accepted
	enddo
4	format(i4, i7, i7, 2f20.13, a15)
40	format(a8, i5, a8, i5, a10, i3, a8, f9.3)

	return
	end subroutine sequence_optimization

	subroutine sequence_optimization_vdwele_dihedral(group, groupdata_backup, Tresista, Tresiend)
	implicit none
	integer*8						:: attempt, i, rotanum_1, rotanum_2, feedback_1
	integer*8						:: ic1, ic2, ip, flag, flag1, stage
	integer*8						:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer*8						:: W_numex(atom_num), W_inb(atom_num,20), W_numex4(atom_num), W_inb4(atom_num, 60)
	integer*8						:: Tresista, Tresiend, Tresinum, rounds
	double precision							:: ran2
	double precision							:: vdw_energy, vdw_energy_min, score, Tscore
	character*4						:: aminoacid_name_1, aminoacid_name_2
	type(groupdetails)				:: group(gnum), Tgroup(gnum)
	type(databackup)				:: groupdata_backup(gnum)
	type(groupdetails)				:: temp_group_1(gnum), temp_group_2(gnum), aa_group_1(40), aa_group_2(40)
	type(energyparameters)			:: tgroup_para(gnum)

	call energy_parameter(group, tgroup_para, gnum) !loads the force field paramters
	call atom_links(group, W_numex, W_inb, W_numex4, W_inb4) !create atom links (again?)

	Tresinum=Tresiend-Tresista+1
	if((Tresinum*3).lt.(sitenum*2)) then
		rounds=Tresinum*3
	else
		rounds=sitenum*2
	endif
	do attempt=1, rounds  !!sitenum is equal to the number amino acids in the peptide that you wanna study
		call ran_gen(ran2,zero)
		ic1=int(ran2*Tresinum-1.0e-3)+1  !! ic1 is randomly picked first
		if(ic1.gt.Tresinum) ic1=Tresinum
		flag = 0
		do while(flag.eq.0)
			call ran_gen(ran2,zero)
			ic2=int(ran2*Tresinum-1.0e-3)+1 !!ic2 is the second amino acid selected for exchange, randomly chosen
			if(ic2.gt.Tresinum) ic2=Tresinum
			if(ic1.ne.ic2) then !!this makes sure that you don't pick the same residue twice
				flag = 1
			endif
		enddo
		ic1=ic1+Tresista-1; ic2=ic2+Tresista-1

		aminoacid_name_1=group(ic1)%gtype
		aminoacid_name_2=group(ic2)%gtype

		call findrotamer(ic1, group, aminoacid_name_1, rotanum_1, aa_group_1)!!calls rotamer possiblities
		call findrotamer(ic2, group, aminoacid_name_2, rotanum_2, aa_group_2)

		flag=0
		if(aminoacid_name_1=="GLY".or.aminoacid_name_1=="ALA".or.aminoacid_name_1=="PRO".or. & 
		aminoacid_name_1=="NGLY".or.aminoacid_name_1=="NALA".or.aminoacid_name_1=="NPRO".or. & 
		aminoacid_name_1=="CGLY".or.aminoacid_name_1=="CALA".or.aminoacid_name_1=="CPRO") then
			vdw_energy_min=1000.0	!initialize
			ip=0 !initialize the rotamer index number
			do i=1, rotanum_1
				call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue in the temporary group

				call check_transplant(zero, ic1, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

				if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
					call vdwenergy_sidechain_receptor(zero, ic1, zero, temp_group_1, tgroup_para, vdw_energy)

					if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
						vdw_energy_min=vdw_energy
						ip=i
					endif
				endif
			enddo

			if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
				call residue_replace(ic1, group, groupdata_backup, ip, aa_group_1, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
				Tgroup=temp_group_1
				group=temp_group_1
				flag=1
			endif
		else

			flag1=0
			do i=1, rotanum_1
				if(i==1) then
					call atom_links4sidechain(ic1, group, S_numex, S_inb, S_numex4, S_inb4)
					temp_group_1=group
					stage=0
					call sidechain_optimization(stage, ic1, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
!					call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
					Tscore=score
					temp_group_2=temp_group_1
				endif

				call residue_replace(ic1, group, groupdata_backup, i, aa_group_1, temp_group_1) !replaces the residue

				stage=0
				call sidechain_optimization(stage, ic1, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
				if(stage==1) then
!					call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
					if(score.le.Tscore) then
						Tscore=score
						temp_group_2=temp_group_1
					endif
					flag1=1
				endif
			enddo

			if(flag1==1) then
				Tgroup=temp_group_2
				group=temp_group_2
				flag=1
			endif
		endif

		if(flag==1) then
			if(aminoacid_name_2=="GLY".or.aminoacid_name_2=="ALA".or.aminoacid_name_2=="PRO".or. & 
			aminoacid_name_2=="NGLY".or.aminoacid_name_2=="NALA".or.aminoacid_name_2=="NPRO".or. & 
			aminoacid_name_2=="CGLY".or.aminoacid_name_2=="CALA".or.aminoacid_name_2=="CPRO") then
				vdw_energy_min=1000.0	!initialize
				ip=0 !initialize the rotamer index number
				do i=1, rotanum_2
					call residue_replace(ic2, Tgroup, groupdata_backup, i, aa_group_2, temp_group_1) !replaces the residue

					call check_transplant(zero, ic2, zero, temp_group_1, feedback_1) !this checks for atomic overlaps and gives answer in feedback_1

					if(feedback_1==1) then !this decision will calculate the VdW energy if there are no atomic overlaps
						call vdwenergy_sidechain_receptor(zero, ic2, zero, temp_group_1, tgroup_para, vdw_energy)

						if(vdw_energy.lt.vdw_energy_min) then !if the calculated energy is lower than the previous minimum, this stores the new minumum and also stores the index to the chosen rotamer
							vdw_energy_min=vdw_energy
							ip=i
						endif
					endif
				enddo

				if(ip>0) then !!this decision checks for an accepted rotamer and effects the actual transplant onto the molecule. then calculates new binding energy
					call residue_replace(ic2, Tgroup, groupdata_backup, ip, aa_group_2, temp_group_1) !group records old PDB file, and temp_group_1 records the new pdb file. the ip tells which rotamer was chosen, and aa_group_1 is all of the rotamer possibilities for the tested/attempted amino acid (ran_resi_1)
					group=temp_group_1
				endif
			else

				flag1=0
				do i=1, rotanum_2
					if(i==1) then
						call atom_links4sidechain(ic2, Tgroup, S_numex, S_inb, S_numex4, S_inb4)
						temp_group_1=Tgroup
						stage=0
						call sidechain_optimization(stage, ic2, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
!						call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
						Tscore=score
						temp_group_2=temp_group_1
					endif

					call residue_replace(ic2, Tgroup, groupdata_backup, i, aa_group_2, temp_group_1) !replaces the residue

					stage=0
					call sidechain_optimization(stage, ic2, temp_group_1, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, score)
					if(stage==1) then
!						call bindingenergy(temp_group_1, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, vdw_new, ele_new, sgb_new, snp_new)
						if(score.le.Tscore) then
							Tscore=score
							temp_group_2=temp_group_1
						endif
						flag1=1
					endif
				enddo

				if(flag1==1) then
					group=temp_group_2
				endif
			endif
		endif
	enddo
	return
	end subroutine sequence_optimization_vdwele_dihedral

	subroutine sequence_optimization_binding(group, Tresista, Tresiend, group_min, & 
		 binding_energy_min, vdw_min, ele_min, sgb_min, snp_min)
	implicit none
	integer*8						:: attempt, feedback_1, feedback_2
	integer*8						:: ic1, ic2
	integer*8						:: W_numex(atom_num), W_inb(atom_num,20), W_numex4(atom_num), W_inb4(atom_num, 60)
	integer*8						:: Tresista, Tresiend, Tresinum, rounds
	double precision							:: ran2
	double precision							:: binding_energy_old, vdw_old, ele_old, sgb_old, snp_old
	double precision							:: binding_energy_new, vdw_new, ele_new, sgb_new, snp_new
	double precision							:: binding_energy_min, vdw_min, ele_min, sgb_min, snp_min
	type(groupdetails)				:: group(gnum), group_min(gnum), group_old(gnum)
	type(groupdetails)				:: Tgroup(gnum), group_new(gnum)
	type(energyparameters)			:: tgroup_para(gnum)

	call energy_parameter(group, tgroup_para, gnum)
	call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
	call bindingenergy(group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)

	group_old=group
	group_min=group_old
	binding_energy_min=binding_energy_old
	vdw_min=vdw_old
	ele_min=ele_old
	sgb_min=sgb_old
	snp_min=snp_old

	Tresinum=Tresiend-Tresista+1
	if((Tresinum*5).lt.(sitenum*4)) then
		rounds=Tresinum*5
	else
		rounds=sitenum*4
	endif
	do attempt=1, rounds  !!sitenum is equal to the number amino acids in the peptide that you wanna study
		call ran_gen(ran2,zero)
		ic1=int(ran2*Tresinum-1.0e-3)+1  !! ic1 is randomly picked first
		if(ic1.gt.Tresinum) ic1=Tresinum
		do while(.true.)
			call ran_gen(ran2,zero)
			ic2=int(ran2*Tresinum-1.0e-3)+1 !!ic2 is the second amino acid selected for exchange, randomly chosen
			if(ic2.gt.Tresinum) ic2=Tresinum
			if(ic1.ne.ic2) then !!this makes sure that you don't pick the same residue twice
				goto 10
			endif
		enddo
10		continue
		ic1=ic1+Tresista-1; ic2=ic2+Tresista-1

		binding_energy_new=binding_energy_old+100.0
		Tgroup=group_old

		feedback_1=0
		if(Tgroup(ic1)%gtype=="GLY".or.Tgroup(ic1)%gtype=="ALA".or.Tgroup(ic1)%gtype=="PRO" & 
		.or.Tgroup(ic1)%gtype=="NGLY".or.Tgroup(ic1)%gtype=="NALA".or.Tgroup(ic1)%gtype=="NPRO".or. & 
		Tgroup(ic1)%gtype=="CGLY".or.Tgroup(ic1)%gtype=="CALA".or.Tgroup(ic1)%gtype=="CPRO") then
			feedback_1=1
		else
			call sidechain_optimization4binding(ic1, Tgroup, tgroup_para, W_numex, W_inb, W_numex4, W_inb4)

			call check_transplant(zero, ic1, zero, Tgroup, feedback_1)
		endif

		feedback_2=0
		if(feedback_1==1) then
			if(Tgroup(ic2)%gtype=="GLY".or.Tgroup(ic2)%gtype=="ALA".or.Tgroup(ic2)%gtype=="PRO" & 
			.or.Tgroup(ic2)%gtype=="NGLY".or.Tgroup(ic2)%gtype=="NALA".or. &
			Tgroup(ic2)%gtype=="NPRO".or.Tgroup(ic2)%gtype=="CGLY".or.Tgroup(ic2)%gtype=="CALA".or. & 
			Tgroup(ic2)%gtype=="CPRO") then
				group_new=Tgroup
				call bindingenergy(group_new, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
								binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
				call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
								ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group_old, group_new, feedback_2)
			else
				group_new=Tgroup
				call sidechain_optimization4binding(ic2, Tgroup, tgroup_para, W_numex, W_inb, W_numex4, W_inb4)
				call check_transplant(zero, ic2, zero, Tgroup, feedback_1)
				if(feedback_1==1) then
					group_new=Tgroup
					call bindingenergy(group_new, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
									binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
					call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
									ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group_old, group_new, feedback_2)
				else
					call bindingenergy(group_new, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, &
									binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
					call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
									ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group_old, group_new, feedback_2)
				endif
			endif
		endif

		if(feedback_2==1) then
			if(binding_energy_old.lt.binding_energy_min) then
				group_min=group_old
				binding_energy_min=binding_energy_old
				vdw_min=vdw_old
				ele_min=ele_old
				sgb_min=sgb_old
				snp_min=snp_old
			endif
		endif

	enddo

	return
	end subroutine sequence_optimization_binding

	subroutine backbone_optimization(group, SA_best_group, groupdata_backup, step, cycle, binding_energy_old, &
		vdw_old, ele_old, sgb_old, snp_old)
	implicit none
	integer*8						:: step, attempt, i, j, k, phipsi_num, cycle
	integer*8						:: feedback_1, feedback_2, flag, flag1, flag2
	integer*8						:: ic(3), ic_min, ic_flag, ran_resi(3)
	integer*8						:: new_min_flag
	double precision							:: binding_energy_old, vdw_old, ele_old, sgb_old, snp_old
	double precision							:: binding_energy_min, vdw_min, ele_min, sgb_min, snp_min
	double precision							:: Tbinding_energy_min, Tvdw_min, Tele_min, Tsgb_min, Tsnp_min
	double precision							:: ran2, com(3)
	type(groupdetails)				:: group(gnum), group_min(gnum), Tgroup_min(gnum), SA_best_group(gnum)
	type(databackup)				:: groupdata_backup(gnum), Tgroupdata_backup(gnum)
	type(databackup)                :: groupdata_backup_min(gnum), groupdata_backup_candidates(20, gnum)
	type(groupdetails)              :: group_candidates(20,gnum)
	type(groupdetails)           	:: Tgroup(gnum)

	new_min_flag = 0

	attempt=1
1	continue

    !decide what kind of CONROT move to do (all internal, one termini, or both termini)
	call ran_gen(ran2,zero)
	if(ran2.le.CONROT_two_term) then
		!two of the pivots for CONROT are the peptide termini
		ic(1)=1
		ic(2)=fragmentnum
		flag = 0
		do while(flag.eq.0)
			call ran_gen(ran2,zero)
			ic(3)=int(ran2*fragmentnum-1.0e-3)+1
			if(ic(3).gt.fragmentnum) ic(3)=fragmentnum
			if(ic(3).ne.ic(1).and.ic(3).ne.ic(2)) then
				flag = 1
			endif
		enddo
	elseif(ran2.le.CONROT_one_term)	then
		!one of the pivots will be a terminus
		call ran_gen(ran2,zero)
		if(ran2.le.0.5) then
			!will use the N-terminus for the terminal pivot
			ic(1)=1
			flag = 0
			do while(flag.eq.0)
				call ran_gen(ran2,zero)
				ic(2)=int(ran2*(fragmentnum-1)-1.0e-3)+2
				if(ic(2).gt.fragmentnum) ic(2)=fragmentnum
				if (ic(2).ne.ic(1)) then
					flag = 1 
				endif
			enddo
			flag = 0
			do while(flag.eq.0)
				call ran_gen(ran2,zero)
				ic(3)=int(ran2*(fragmentnum-1)-1.0e-3)+2
				if(ic(3).gt.fragmentnum) ic(3)=fragmentnum
				if(ic(3).ne.ic(1).and.ic(3).ne.ic(2)) then
					flag = 1
				endif
			enddo
		else
			!will use the C-terminus for the terminal pivot
			ic(1)=fragmentnum
			flag = 0
			do while(flag.eq.0)
				call ran_gen(ran2,zero)
				ic(2)=int(ran2*(fragmentnum-1)-1.0e-3)+1
				if(ic(2).eq.fragmentnum) ic(2)=fragmentnum-1
				if (ic(2).ne.ic(1)) then
					flag=1
				endif
			enddo
			flag = 0
			do while(flag.eq.0)
				call ran_gen(ran2,zero)
				ic(3)=int(ran2*(fragmentnum-1)-1.0e-3)+1
				if(ic(3).eq.fragmentnum) ic(3)=fragmentnum-1
				if(ic(3).ne.ic(1).and.ic(3).ne.ic(2)) then
					flag = 1
				endif
			enddo
		endif

	else
		!pick three internal peptide residues randomly
		call ran_gen(ran2,zero)
		ic(1)=int(ran2*(fragmentnum-2)-1.0e-3)+2
		if(ic(1).eq.fragmentnum) ic(1)=fragmentnum-1
		flag = 0
		do while(flag.eq.0)
			call ran_gen(ran2,zero)
			ic(2)=int(ran2*(fragmentnum-2)-1.0e-3)+2
			if(ic(2).eq.fragmentnum) ic(2)=fragmentnum-1
			if (ic(2).ne.ic(1)) then
				flag = 1
			endif
		enddo
		flag = 0
		do while(flag.eq.0)
			call ran_gen(ran2,zero)
			ic(3)= int(ran2*(fragmentnum-2)-1.0e-3)+2
			if(ic(3).eq.fragmentnum) ic(3)=fragmentnum-1
			if(ic(3).ne.ic(1).and.ic(3).ne.ic(2)) then
				flag=1
			endif
		enddo
	endif

    !order the three residues in ic1
	do i=1, 2
		ic_min=ic(i)
		ic_flag=i
		do j=i+1, 3
			if(ic_min.gt.ic(j)) then
				ic_min=ic(j)
				ic_flag=j
			endif
		enddo
		ic(ic_flag)=ic(i)
		ic(i)=ic_min
	enddo

	do i=1, 3
		ran_resi(i)=residuesite(ic(i))
	enddo

	flag1=0
	flag2=0
	if(ran_resi(1)==1) flag1=1
	if(ran_resi(3)==sitenum) flag2=1

	if(flag1.ne.1.and.flag2.ne.1) then
		CONF_Internal_attempts = CONF_Internal_attempts+1
		call concertedrotation_center(group, groupdata_backup, ran_resi, phipsi_num, group_candidates, groupdata_backup_candidates)

		if(phipsi_num.eq.0) then
			open(3, file="energydetails.txt", access="append")
				write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "CRA-mid"
				!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
				write(3,*) (group(j)%gtype, j=1, sitenum)
				write(3,*) "*******************************"
			close(3)

		else
			binding_energy_min=100.0
			do k=1, phipsi_num
				do j=1, gnum
					Tgroup(j)=group_candidates(k,j)
				enddo
				do j=1, gnum
					Tgroupdata_backup(j)%coo=groupdata_backup_candidates(k,j)%coo
				enddo

				call check_backbone(Tgroup, feedback_1)

				!check backbone change doesn't make peptide bend too far from surface
				com = calc_COM(group, 1)

				if(feedback_1.ne.0 .and. (com(3)-start_pep_com(3)).lt.max_zdist_from_surface) then
					call sequence_optimization_vdwele_dihedral(Tgroup, Tgroupdata_backup, ran_resi(1), ran_resi(3))
					call check_transplant(one+1, zero, zero, Tgroup, feedback_1)
					if(feedback_1.ne.0) then
						call sequence_optimization_binding(Tgroup, ran_resi(1), ran_resi(3), &
														Tgroup_min, Tbinding_energy_min, Tvdw_min, &
														Tele_min, Tsgb_min, Tsnp_min)
						if(Tbinding_energy_min.lt.binding_energy_min) then
							group_min=Tgroup_min
							binding_energy_min=Tbinding_energy_min
							vdw_min=Tvdw_min
							ele_min=Tele_min
							sgb_min=Tsgb_min
							snp_min=Tsnp_min
							groupdata_backup_min=Tgroupdata_backup
						endif
					endif
				endif

			enddo

			call MC_technique_backbone(binding_energy_min, binding_energy_old, vdw_min, vdw_old, ele_min, ele_old, &
				sgb_min, sgb_old, snp_min, snp_old, group, group_min, feedback_2)
			if(feedback_2==1) then
				CONF_Internal_success = CONF_Internal_success+1
				flag4conformer=1
				groupdata_backup=groupdata_backup_min
				rmsd_flag = 1

				call generatepdb(cycle, step, attempt, group)
				
				if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
					call Update_Top_Peptides(binding_energy_old, cycle, step, one)
				endif

				if (binding_energy_old.lt.energy_min-.01) then
					energy_min = binding_energy_old
					if (SA_flag.eq.1) then
						SA_best_group = group
						SA_best_hyd = cur_hyd
						SA_above_hyd = above_hyd
						SA_below_hyd = below_hyd
						SA_trp_count = trp_count
					endif
					new_min_flag = 1
				endif
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "CRA-mid"
				!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
				write(3,*) (group(j)%gtype, j=1, sitenum)
				write(3,*) "*******************************"
			close(3)
		endif

	elseif(flag1.eq.1.and.flag2.eq.1) then
		CONF_N2C_attempts = CONF_N2C_attempts + 1
		call concertedrotation_whole(group, groupdata_backup, ran_resi, phipsi_num, group_candidates, groupdata_backup_candidates)

		if(phipsi_num.eq.0) then
			goto 1
		endif

		binding_energy_min=100.0
		do k=1, phipsi_num
			do j=1, gnum
				Tgroup(j)=group_candidates(k,j)
			enddo
			do j=1, gnum
				Tgroupdata_backup(j)%coo=groupdata_backup_candidates(k,j)%coo
			enddo

			call check_backbone(Tgroup, feedback_1)

			!check backbone change doesn't make peptide bend too far from surface
			com = calc_COM(group, 1)

			if(feedback_1.eq.0 .or. (com(3)-start_pep_com(3)).ge.max_zdist_from_surface)  goto 40

			call sequence_optimization_vdwele_dihedral(Tgroup, Tgroupdata_backup, ran_resi(1), ran_resi(3))

			call check_transplant(one+1, zero, zero, Tgroup, feedback_1)

			if(feedback_1==0) goto 40

			call sequence_optimization_binding(Tgroup, ran_resi(1), ran_resi(3), &
												Tgroup_min, Tbinding_energy_min, Tvdw_min, Tele_min, Tsgb_min, Tsnp_min)

			if(Tbinding_energy_min.lt.binding_energy_min) then
				group_min=Tgroup_min
				binding_energy_min=Tbinding_energy_min
				vdw_min=Tvdw_min
				ele_min=Tele_min
				sgb_min=Tsgb_min
				snp_min=Tsnp_min
				groupdata_backup_min=Tgroupdata_backup
			endif
40			continue
		enddo

		call MC_technique_backbone(binding_energy_min, binding_energy_old, vdw_min, vdw_old, ele_min, ele_old, &
			 sgb_min, sgb_old, snp_min, snp_old, group, group_min, feedback_2)
		if(feedback_2==1) then
			CONF_N2C_success = CONF_N2C_success + 1
			flag4conformer=1
			rmsd_flag = 1
			groupdata_backup=groupdata_backup_min

			if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
				call Update_Top_Peptides(binding_energy_old, cycle, step, one)
			endif

			if (binding_energy_old.lt.energy_min-.01) then
				energy_min = binding_energy_old
				if (SA_flag.eq.1) then
					SA_best_group = group
					SA_best_hyd = cur_hyd
					SA_above_hyd = above_hyd
					SA_below_hyd = below_hyd
					SA_trp_count = trp_count
				endif
				new_min_flag = 1
			endif

			call generatepdb(cycle, step, attempt, group)
		endif

		open(3, file="energydetails.txt",  access="append")
			write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "CRA-full"
			!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
			write(3,*) (group(j)%gtype, j=1, sitenum)
			write(3,*) "*******************************"
		close(3)

	elseif(flag1.eq.1.and.flag2.ne.1) then
		CONF_N_attempts = CONF_N_attempts + 1
		call concertedrotation_Nterm(group, groupdata_backup, ran_resi, phipsi_num, group_candidates, groupdata_backup_candidates)

		if(phipsi_num.eq.0) then
			goto 1
		endif

		binding_energy_min=100.0
		do k=1, phipsi_num
			do j=1, gnum
				Tgroup(j)=group_candidates(k,j)
			enddo
			do j=1, gnum
				Tgroupdata_backup(j)%coo=groupdata_backup_candidates(k,j)%coo
			enddo

			call check_backbone(Tgroup, feedback_1)

			!check backbone change doesn't make peptide bend too far from surface
			com = calc_COM(group, 1)

			if(feedback_1.eq.0 .or. (com(3)-start_pep_com(3)).ge.max_zdist_from_surface) goto 50

			call sequence_optimization_vdwele_dihedral(Tgroup, Tgroupdata_backup, ran_resi(1), ran_resi(3))

			call check_transplant(one+1, zero, zero, Tgroup, feedback_1)

			if(feedback_1==0) goto 50

			call sequence_optimization_binding(Tgroup, ran_resi(1), ran_resi(3), &
												Tgroup_min, Tbinding_energy_min, Tvdw_min, Tele_min, Tsgb_min, Tsnp_min)

			if(Tbinding_energy_min.lt.binding_energy_min) then
				group_min=Tgroup_min
				binding_energy_min=Tbinding_energy_min
				vdw_min=Tvdw_min
				ele_min=Tele_min
				sgb_min=Tsgb_min
				snp_min=Tsnp_min
				groupdata_backup_min=Tgroupdata_backup
			endif
50			continue
		enddo

		call MC_technique_backbone(binding_energy_min, binding_energy_old, vdw_min, vdw_old, ele_min, ele_old, &
			 sgb_min, sgb_old, snp_min, snp_old, group, group_min, feedback_2)
		if(feedback_2==1) then
			CONF_N_success = CONF_N_success + 1
			flag4conformer=1
			groupdata_backup=groupdata_backup_min
			rmsd_flag = 1
			call generatepdb(cycle, step, attempt, group)

			if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
				call Update_Top_Peptides(binding_energy_old, cycle, step, one)
			endif

			if (binding_energy_old.lt.energy_min-.01) then
				energy_min = binding_energy_old
				if (SA_flag.eq.1) then
					SA_best_group = group
					SA_best_hyd = cur_hyd
					SA_above_hyd = above_hyd
					SA_below_hyd = below_hyd
					SA_trp_count = trp_count
				endif
				new_min_flag = 1
			endif
		endif

		open(3, file="energydetails.txt",  access="append")
			write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "CRA-Nend"
			!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
			write(3,*) (group(j)%gtype, j=1, sitenum)
			write(3,*) "*******************************"
		close(3)

	elseif(flag1.ne.1.and.flag2.eq.1) then
		CONF_C_attempts = CONF_C_attempts + 1
		call concertedrotation_Cterm(group, groupdata_backup, ran_resi, phipsi_num, group_candidates, groupdata_backup_candidates)

		if(phipsi_num.eq.0) then
			goto 1
		endif

		binding_energy_min=100.0
		do k=1, phipsi_num
			do j=1, gnum
				Tgroup(j)=group_candidates(k,j)
			enddo
			do j=1, gnum
				Tgroupdata_backup(j)%coo=groupdata_backup_candidates(k,j)%coo
			enddo

			call check_backbone(Tgroup, feedback_1)

			!check backbone change doesn't make peptide bend too far from surface
			com = calc_COM(group, 1)

			if(feedback_1.eq.0 .or. (com(3)-start_pep_com(3)).ge.max_zdist_from_surface)  goto 60

			call sequence_optimization_vdwele_dihedral(Tgroup, Tgroupdata_backup, ran_resi(1), ran_resi(3))

			call check_transplant(one+1, zero, zero, Tgroup, feedback_1)

			if(feedback_1==0) goto 60

			call sequence_optimization_binding(Tgroup, ran_resi(1), ran_resi(3), &
												Tgroup_min, Tbinding_energy_min, Tvdw_min, Tele_min, Tsgb_min, Tsnp_min)

			if(Tbinding_energy_min.lt.binding_energy_min) then
				group_min=Tgroup_min
				binding_energy_min=Tbinding_energy_min
				vdw_min=Tvdw_min
				ele_min=Tele_min
				sgb_min=Tsgb_min
				snp_min=Tsnp_min
				groupdata_backup_min=Tgroupdata_backup
			endif
60			continue
		enddo

		call MC_technique_backbone(binding_energy_min, binding_energy_old, vdw_min, vdw_old, ele_min, ele_old, &
			 sgb_min, sgb_old, snp_min, snp_old, group, group_min, feedback_2)
		if(feedback_2==1) then
			CONF_C_success = CONF_C_success + 1
			flag4conformer=1
			groupdata_backup=groupdata_backup_min
			rmsd_flag=1
			call generatepdb(cycle, step, attempt, group)

			if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
				call Update_Top_Peptides(binding_energy_old, cycle, step, one)
			endif

			if (binding_energy_old.lt.energy_min-.01) then
				energy_min = binding_energy_old
				if (SA_flag.eq.1) then
					 SA_best_group = group
					 SA_best_hyd = cur_hyd
					 SA_above_hyd = above_hyd
					 SA_below_hyd = below_hyd
					 SA_trp_count = trp_count
				endif
				new_min_flag = 1
			endif
		endif

		open(3, file="energydetails.txt",  access="append")
			write(3,4) cycle, step, attempt, binding_energy_old, (vdw_old+ele_old+sgb_old), "CRA-Cend"
			!write(3,"(<sitenum>(a4))") (group(j)%gtype, j=1, sitenum)
			write(3,*) (group(j)%gtype, j=1, sitenum)
			write(3,*) "*******************************"
		close(3)

	endif

	if(new_min_flag.eq.1) then !if we found a new minimum, then store the energy data
		open(2, file="minimum_energy.txt", access="append")
			write(2,41) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
		close(2)
	endif
5	continue
4	format(i5,i7,i7,2f20.13,a15)
41	format(a8, i5, a8, i5, a10, i3, a8, f9.3)

	return
	end subroutine backbone_optimization

	subroutine calc_mutate_ps(energy_per_res, res_select_p_vector)
	implicit none
	double precision		:: energy_per_res(sitenum,4)
	double precision		:: res_select_p_vector(sitenum), net_P, cur_P
	integer*8	:: i

	!will use Boltzmann weighting to determine which residue mutation probabilities
	net_P = 0
	do i=1, sitenum
		if (decomp_score_type.eq."REG") then
			cur_P = exp(-energy_per_res(i,4)/ekt_seq)
		else
			cur_P = exp(energy_per_res(i,4)/ekt_seq)
		endif
		res_select_p_vector(i) = cur_P + net_P
		net_P = net_P + cur_P
	enddo
	res_select_p_vector = res_select_p_vector/net_P

	end subroutine calc_mutate_ps

	subroutine update_TABU_ban_list()
	implicit none
	
	integer	::	i, TABU_swap

	do i=1,sitenum
		TABU_ban_lengths(i) = TABU_ban_lengths(i) - 1
		if (TABU_ban_lengths(i).eq.0) then
			TABU_allowed_length = TABU_allowed_length + 1

			TABU_swap = TABU_allowed_residues(TABU_allowed_length)
			TABU_allowed_residues(TABU_allowed_length) = TABU_allowed_residues(i)
			TABU_allowed_residues(i) = TABU_swap

			TABU_swap = TABU_ban_lengths(TABU_allowed_length)
			TABU_ban_lengths(TABU_allowed_length) = TABU_ban_lengths(i)
			TABU_ban_lengths(i) = TABU_swap
		endif
	enddo
	end subroutine update_TABU_ban_list

	subroutine Update_Top_Peptides(cur_score, cycle, step, attempt)
	implicit none
	double precision			:: cur_score, worst_score
	character*30	:: descriptor
	character*5		:: cyclechar, stepchar, attemptchar
	integer*8		:: worst_index, cycle, step,attempt

	write(cyclechar,"(i3)") cycle
	write(stepchar,"(i5)") step
	write(attemptchar,"(i4)") attempt
	descriptor = trim(adjustl(cyclechar))//'_'//trim(adjustl(stepchar))//'_'//trim(adjustl(attemptchar))//".pdb"

	!update score array
	top_pep_scores(ESURF_worst_index) = cur_score
	ESURF_pep_descriptors(ESURF_worst_index) = descriptor

	!find worst score in new array
	call Find_Worst_Index(one, ESURF_count, worst_score, worst_index, one)
	ESURF_worst_score = worst_score
	ESURF_worst_index = worst_index

	end subroutine Update_Top_Peptides

	subroutine Calc_Score_Smart_Restart(smart_restart_scores, smart_restart_best_score, smart_restart_cur_score, smart_restart_cur_flag, smart_restart_cur_len)
	implicit none
	double precision          :: smart_restart_scores(smart_restart_len), smart_restart_best_score, smart_restart_cur_score
	integer*8		:: smart_restart_cur_flag, smart_restart_cur_len

	if (smart_restart_score_method.eq.1.or.smart_restart_score_method.eq.3) then
		!calculate average score of current interval
		smart_restart_cur_score = sum(smart_restart_scores)/smart_restart_cur_len !NOTE: need to use "cur_len" variable rather than "len" variable for properly calculating score for interval without initial randomization. Either of variables will work for other cycles
	elseif (smart_restart_score_method.eq.2) then
		!calculate best score in current interval
		smart_restart_cur_score = minval(smart_restart_scores(1:smart_restart_cur_len)) !same comment as above about "cur_len" vs. "len" variables
	else 
		open(20, file="error.txt", access="append")
        write(20,*) "Unable to use current smart_restart_score_method of ", smart_restart_score_method
		write(20,*) "This should never happen as the error should have been caught when reading the input file."
		write(20,*) "Terminating program"
        close(20)
        stop
	endif

	if (smart_restart_cur_score.gt.(smart_restart_best_score+smart_restart_offset)) then
		smart_restart_cur_flag = 1 !current cycle not doing so great, so do a random restart
	elseif (smart_restart_cur_score.lt.smart_restart_best_score) then
		!current cycle doing better than previous best cycle, so update best score
		if (smart_restart_score_method.eq.3) then
			smart_restart_best_score = minval(smart_restart_scores)
		else
			smart_restart_best_score = smart_restart_cur_score
		endif
	endif

	end subroutine Calc_Score_Smart_Restart

	subroutine Smart_Restart_Update(smart_restart_cur_score, smart_restart_scores, smart_restart_index, binding_energy_old, smart_restart_cur_len, &
									smart_restart_cur_flag, smart_restart_best_score, smart_restart_full_flag)
	implicit none
	double precision          :: smart_restart_scores(smart_restart_len)
	double precision			:: binding_energy_old, smart_restart_cur_score, smart_restart_best_score
	integer*8		:: smart_restart_index, smart_restart_full_flag, smart_restart_cur_len, smart_restart_cur_flag

		smart_restart_scores(smart_restart_index) = binding_energy_old
		smart_restart_index = mod(smart_restart_index, smart_restart_len) + 1
		if (smart_restart_full_flag.eq.0) then
			!smart restart score array not yet filled up, so add value and keep going
			smart_restart_cur_len = smart_restart_cur_len + 1
			if (smart_restart_cur_len.eq.smart_restart_len) then
				smart_restart_full_flag = 1
				call Calc_Score_Smart_Restart(smart_restart_scores, smart_restart_best_score, smart_restart_cur_score, smart_restart_cur_flag,  smart_restart_cur_len)
			endif
		else
			!We have filled up the array of scores, so average the values and check if we should restart design
			call Calc_Score_Smart_Restart(smart_restart_scores, smart_restart_best_score, smart_restart_cur_score, smart_restart_cur_flag, smart_restart_cur_len)
		endif

	end subroutine Smart_Restart_Update

    subroutine Get_Min_Max_Box(i, res_min, res_max, group)
    implicit none
    integer *8      :: i, j
    double precision      :: res_min(3), res_max(3), cur_coords(3)
    type(groupdetails)				:: group(gnum)

    res_min = 10000
    res_max = -10000
    do j=1,group(1)%cnum1
        cur_coords = group(i)%coo1(j,:)
        if (cur_coords(1).lt.res_min(1)) res_min(1) = cur_coords(1)
        if (cur_coords(2).lt.res_min(2)) res_min(2) = cur_coords(2)
        if (cur_coords(3).lt.res_min(3)) res_min(3) = cur_coords(3)
        if (cur_coords(1).gt.res_max(1)) res_max(1) = cur_coords(1)
        if (cur_coords(2).gt.res_max(2)) res_max(2) = cur_coords(2)
        if (cur_coords(3).gt.res_max(3)) res_max(3) = cur_coords(3)
    enddo

    do j=1,group(1)%cnum2
        cur_coords = group(i)%coo2(j,:)
        if (cur_coords(1).lt.res_min(1)) res_min(1) = cur_coords(1)
        if (cur_coords(2).lt.res_min(2)) res_min(2) = cur_coords(2)
        if (cur_coords(3).lt.res_min(3)) res_min(3) = cur_coords(3)
        if (cur_coords(1).gt.res_max(1)) res_max(1) = cur_coords(1)
        if (cur_coords(2).gt.res_max(2)) res_max(2) = cur_coords(2)
        if (cur_coords(3).gt.res_max(3)) res_max(3) = cur_coords(3)
    enddo

    do j=1,group(1)%cnum3
        cur_coords = group(i)%coo3(j,:)
        if (cur_coords(1).lt.res_min(1)) res_min(1) = cur_coords(1)
        if (cur_coords(2).lt.res_min(2)) res_min(2) = cur_coords(2)
        if (cur_coords(3).lt.res_min(3)) res_min(3) = cur_coords(3)
        if (cur_coords(1).gt.res_max(1)) res_max(1) = cur_coords(1)
        if (cur_coords(2).gt.res_max(2)) res_max(2) = cur_coords(2)
        if (cur_coords(3).gt.res_max(3)) res_max(3) = cur_coords(3)
    enddo

    end subroutine Get_Min_Max_Box

    subroutine find_Neighbors(start, stop, group)
    implicit none

    integer*8   :: i, j, index, start, stop, pol_i
    double precision      :: res_min(3), res_max(3), x_dis, y_dis, z_dis, net_dist
    type(groupdetails)				:: group(gnum)

    !this loops finds the neighbors for all ligand residues
    do i=start,stop
        index = 1
        !get min/max coordinates of the current residue
        call Get_Min_Max_Box(index, res_min, res_max, group)
        !iterate through all polymer residues to see which are within cutoff 
        do j=sitenum+1, gnum
            pol_i = j-sitenum !shift index to go to correct entry in neighbor_rec_res_box array (this array starts at 1 for the first receptor residue)
            x_dis = min(abs(res_min(1)-neighbor_rec_res_box_edges(pol_i,1)), abs(res_min(1)-neighbor_rec_res_box_edges(pol_i,4)), &
                        abs(res_max(1)-neighbor_rec_res_box_edges(pol_i,1)), abs(res_max(1)-neighbor_rec_res_box_edges(pol_i,4)))
            if (x_dis.lt.neighbor_cutoff) then
                y_dis = min(abs(res_min(2)-neighbor_rec_res_box_edges(pol_i,2)), abs(res_min(2)-neighbor_rec_res_box_edges(pol_i,5)), &
                        abs(res_max(2)-neighbor_rec_res_box_edges(pol_i,2)), abs(res_max(2)-neighbor_rec_res_box_edges(pol_i,5)))
                net_dist = x_dis * x_dis + y_dis * y_dis
                if (net_dist.lt.neighbor_cutoff2) then
                    z_dis = min(abs(res_min(3)-neighbor_rec_res_box_edges(pol_i,3)), abs(res_min(3)-neighbor_rec_res_box_edges(pol_i,6)), &
                    abs(res_max(3)-neighbor_rec_res_box_edges(pol_i,3)), abs(res_max(3)-neighbor_rec_res_box_edges(pol_i,6)))
                    net_dist = net_dist + z_dis * z_dis
                    if (net_dist.lt.neighbor_cutoff2) then
                        neighbor_list(i, index) = j
                        index = index + 1
                    endif
                endif
            endif
        enddo
        neighbor_list_lens(i) = index - 1
    enddo

    end subroutine find_Neighbors

    subroutine Get_Receptor_Res_Boxes(group)
    implicit none

    type(groupdetails)				:: group(gnum)
    integer*8   :: i, index
    double precision      :: res_min(3), res_max(3)

    index = 1
    do i=sitenum+1, gnum
        call Get_Min_Max_Box(index, res_min, res_max, group)
        neighbor_rec_res_box_edges(i,1:3) = res_min
        neighbor_rec_res_box_edges(i,4:6) = res_min
    enddo
    
    end subroutine Get_Receptor_Res_Boxes

    subroutine print_Neighbors()
    implicit none
    integer*8          :: i,j

    open(10, file="Neighbors.txt")
    do i=1,sitenum
        write(10,*) "Neighbors of peptide residue", i
        do j=1, neighbor_list_lens(i)
            write(10,*) neighbor_list(i,j)
        enddo
        write(10,*) ""; write(10,*), ""
    enddo
    close(10)

    end subroutine print_Neighbors

	subroutine rigid_body_move(group, SA_best_group, groupdata_backup, step, cycle, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old, group_para)
	implicit none
	
	integer*8					:: step, attempt, j, cycle
	integer*8					:: feedback
	integer*8					:: move_type
	double precision						:: binding_energy_old, vdw_old, ele_old, sgb_old, snp_old
	double precision						:: binding_energy_new, vdw_new, ele_new, sgb_new, snp_new
	double precision						:: ran2
	integer*8 					:: W_numex(atom_num), W_numex4(atom_num)
	integer*8 					:: W_inb(atom_num,20), W_inb4(atom_num, 60)
	type(groupdetails)			:: group(gnum), Tgroup(gnum), SA_best_group(gnum)
	type(databackup)			:: groupdata_backup(gnum) 
	type(energyparameters) 		:: group_para(gnum)

	Tgroup = group
	
	!pick either a translation or rotation rigid body move
	call ran_gen(ran2,zero)
	if (ran2.lt.rotation_switch) then
		call random_ligand_translation(Tgroup)
		move_type=1
	else
		call random_ligand_rotation(Tgroup)
		move_type=2
	endif

	!create atom links for proper evaluation of binding energy
	call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)

	!calculate new binding energy and determine if rigid body move should be accepted
	call bindingenergy(Tgroup, group_para, W_numex, W_inb, W_numex4, W_inb4, &
	binding_energy_new, vdw_new, ele_new, sgb_new, snp_new)
	call MC_technique(binding_energy_new, binding_energy_old, vdw_new, vdw_old, &
   	ele_new, ele_old, sgb_new, sgb_old, snp_new, snp_old, group, Tgroup, feedback)

	if (feedback.eq.1) then !change accepted
		!update peptide center of mass position (only necessary for translation, since rotations are done about peptide center of mass)
		!cur_pep_com = cur_pep_com + translate_vec !currently not using this option, since center of mass not updated after backbone or sequence changes

		!update backup coordinates (necessary for mutations involving proline)
		call update_backup_coords(group, groupdata_backup)

		!write a pdb file
		call generatepdb(cycle, step, one, group)

		!check if peptide should be added to final ESURF calculation
		if (ESURF_flag.eq.1.and.binding_energy_old.lt.ESURF_worst_score) then
			call Update_Top_Peptides(binding_energy_old, cycle, step, one)
		endif

		!check if we found the minimum energy during design
		if (binding_energy_old.lt.energy_min-.01) then
			energy_min = binding_energy_old
			if (SA_flag.eq.1) then
				SA_best_group = group
				SA_best_hyd = cur_hyd
				SA_above_hyd = above_hyd
				SA_below_hyd = below_hyd
				SA_trp_count = trp_count
			endif
			open(2, file="minimum_energy.txt", access="append")
				write(2,40) "cycle", cycle, "step", step, "attempt", attempt, "EBind", binding_energy_old
			close(2)
		endif
	endif	
40  format(a8, i5, a8, i5, a10, i3, a8, f9.3)

	open(3, file="energydetails.txt",  access="append")
		if(move_type.eq.1) then
			write(3,4) cycle, step, binding_energy_old, (vdw_old+ele_old+sgb_old), "        TRANSLATE" 
		else
			write(3,4) cycle, step, binding_energy_old, (vdw_old+ele_old+sgb_old), "        ROTATE" 
		endif
		write(3,*) (group(j)%gtype, j=1, sitenum)
		write(3,*) "*******************************"
	close(3)
4	format(i4, i7, 2f20.13, a15)

	end subroutine rigid_body_move

end module optimization_techniques

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program ProteinDesign

	use constant
	use randomgenerator
	use input
	use pdbfile
	use mathfunction
	use database
	use energy_calculation
	use advancedfunction
	use optimization_techniques
    use omp_lib

	implicit none
	integer*8						:: step, i, sub_circle, stat
	double precision							:: ran2, rmsd=0.0
	double precision							:: binding_energy_old, vdw_old, ele_old, sgb_old, snp_old
    double precision					        :: smart_restart_cur_score=1000.0, smart_restart_best_score=1000.0  !stores score over previous len steps for current cycle, and current best score
    double precision, allocatable             :: smart_restart_scores(:)
    integer*8                       :: smart_restart_index = 1, smart_restart_cur_len = 0, smart_restart_full_flag = 0, smart_restart_cur_flag = 0
	integer*8						:: SA_cycle = 1, net_step = 1
	integer*8 						:: W_numex(atom_num), W_numex4(atom_num)
	integer*8 						:: W_inb(atom_num,20), W_inb4(atom_num, 60)
	double precision							:: volume, surfarea, comp_esurf, rec_esurf, lig_esurf, net_esurf
	character*30					:: fname !used for loading pdbfiles during ESURF calculations
	type(databackup)				:: groupdata_backup(gnum) 
	type(groupdetails)				:: group(gnum)
	type(groupdetails)				:: SA_best_group(gnum) !stores best peptide found during a simulated annealing interval, which is passed along to next round
	type(groupdetails) 				:: temp_group(gnum)
	type(energyparameters) 			:: group_para(gnum)

	open(5, file="energydetails.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="minimum_energy.txt",status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="randomnumber.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="rmsd.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

    open(5, file="error.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="output.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="energydecomp.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="ESURF_Results.txt", status="old", iostat=stat)
	if (stat == 0) close(5, status="delete")

	open(5, file="energyprofile.txt", status="replace")
			write(5,*) "Subcycle   Step   binding_Energy   VDW_Energy   ELE+EGB_Energy  Scaled_Peptide_Energy"
	close(5)

	open(5, file="energydecomp.txt", status="replace")
		write(5,*) "Resname   VDW    ELE    EGB    NET"
	close(5)
	
	call inputfile
    !call allocate_memory(filename1, group, groupdata_backup, gnum, atom_num, W_numex, W_numex4, W_inb, W_inb4, group_para)

	allocate(top_pep_scores(num_top_layer)) !currently reusing the top_pep_scores variable to calculate top of polymer z surface, since it is needed by Find_Worst_Index function
										  ! Doesn't cause problems, but not the nicest looking solution. Perhaps prettify this in future
	top_pep_scores = -1000
	call readpdb(group, groupdata_backup, gnum)
	deallocate(top_pep_scores)

	call PH_checking(group, gnum)
	call rotamerlib
	call ramachandranmap

	!move peptide over top of surface, if applicable
	if (move_pep_over_surface_flag.eq.1) then
		call move_pep_over_surface(group)
	else
		start_pep_com = calc_COM(group, 1) !note that move_pep_over_surface already calculates the center of mass, so don't need to do it again
	endif

	allocate(energy_per_res(sitenum, 4))
	!allocate(res_select_p_vector(sitenum))
	allocate(SEQ_PerRes_attempts(sitenum))
	allocate(SEQ_PerRes_success(sitenum))
	SEQ_PerRes_attempts = 0; SEQ_PerRes_success = 0
    !if (neighbor_flag.eq.1) then
        !*** below variables are currently not used! ***
        !allocate(neighbor_list(sitenum, neighbor_length))
        !allocate(neighbor_list_lens(sitenum))
        !allocate(neighbor_rec_res_box_edges(gnum-sitenum,6))

        !get boxes for receptor residues
        !call Get_Receptor_Res_Boxes(group)
        !call print_Neighbors()
        !stop
    !endif 

	!set up array for storing top peptides to calculate ESURF on, if applicable
	if (ESURF_flag.eq.1) then
		allocate(top_pep_scores(ESURF_count))
		top_pep_scores = 1000 !initialze score array to very high number
		ESURF_worst_index = 1
		allocate(ESURF_pep_descriptors(ESURF_count))
	endif

	!set up likelihood for mutating residues
	!if (decomp_mutate_flag.eq.1) then
    !do i=1, sitenum
    !    res_select_p_vector(i) = 1.0*i/sitenum
    !enddo
	!endif

	!setup array for storing scores in smart restart, if applicable
    !if (smart_restart_flag.eq.1) then
    !    allocate(smart_restart_scores(smart_restart_len))
    !endif

    if(recalcu_switch==0) then
        call scmf_substitution(group, groupdata_backup, sub_circle, temp_group) !this adjusts hydration properties of peptide to match input file specifications

		!replace tryptophans with valine if above the tryptophan limit
		call replace_extra_tryptophans(group, groupdata_backup, temp_group)
	endif

	rec_calculation_flag=1 !this label says that no calculation has been done for the receptor
	if (SA_flag.eq.1) SA_best_group = group !if doing simulated annealing, then store best group for current round as the current group

	!*********************************
	! THE INITIAL RANDOMIZATION STAGE
	!*********************************

	if (initial_randomize.eq.0) then
		!do a short design before randomizing peptide sequence to get baseline for a "good" peptide binder
		call energy_parameter(group, group_para, gnum) !this loads the amber force field parameters for the pdb info into the program
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
		call bindingenergy(group, group_para, W_numex, W_inb, W_numex4, W_inb4, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
		energy_min = binding_energy_old 
		open(5, file="energyprofile.txt", access="append") !creates output file for energy calculations
			write(5,6) 0, 0, binding_energy_old, vdw_old, (ele_old+sgb_old), snp_old
		close(5)

		!NOTE if are doing simulated annealing, then we use the final temperature for this step
		if (SA_flag.eq.1) then
			ekt_seq = SA_T_end
			ekt_backbone = SA_T_end
		endif

		step=1
		do while (step.le.steps_before_randomize)
            call ran_gen(ran2,zero)
            if(ran2.ge.seq_switch) then
                !change peptide sequence
                call sequence_optimization(group, SA_best_group, groupdata_backup, step, zero, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
            elseif (ran2.ge.backbone_switch) then
                !change peptide backbone
                call backbone_optimization(group, SA_best_group, groupdata_backup, step, zero, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
			else
				!rotate or translate peptide
				call rigid_body_move(group, SA_best_group, groupdata_backup, step, zero, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old, group_para)
			endif

			!update mutation probability array, if applicable
			!if (decomp_mutate_flag.eq.1) then
			!	call calc_mutate_ps(energy_per_res, res_select_p_vector)
			!endif

			open(5, file="energyprofile.txt", access="append")
				write(5,6) 0, step, binding_energy_old, vdw_old, (ele_old+sgb_old), snp_old
			close(5)

			open(5, file="energydecomp.txt", access="append")
			write(5,702) "cycle", zero, "step", step 
			do i=1,sitenum
				write(5, 701) group(i)%gtype, energy_per_res(i,1), energy_per_res(i,2), energy_per_res(i,3), energy_per_res(i,4)
			enddo

			open(5, file="rmsd.txt", access="append")
				if (rmsd_flag.eq.1) then
					call rmsd_calculation(group, rmsd)
					rmsd_flag = 0
					write(5,703) zero, step, rmsd
				endif
				write(5,*) step, rmsd
			close(5)

			!update dynamic cycle info
            !if (smart_restart_flag.eq.1) then
			!	call Smart_Restart_Update(smart_restart_cur_score, smart_restart_scores, smart_restart_index, binding_energy_old, &
			!	smart_restart_cur_len, smart_restart_cur_flag, smart_restart_best_score, smart_restart_full_flag)
			!	if (smart_restart_cur_flag.eq.1) then
			!		!current cycle not doing so hot, so do a random restart
			!		step = steps_before_randomize
			!	endif
            !endif
			step = step+1
		close(5)
		enddo

		!if (smart_restart_flag.eq.1) then
		!	!reset smart restart values before beginning design
		!	if (smart_restart_best_score.eq.1000.0) then
		!		!never calculated score during intial design, so calculate it now
		!		call Smart_Restart_Update(smart_restart_cur_score, smart_restart_scores, smart_restart_index, binding_energy_old, &
		!		smart_restart_cur_len, smart_restart_cur_flag, smart_restart_best_score, one)
		!	endif
		!	smart_restart_cur_len = 0
		!	smart_restart_index = 1
		!	smart_restart_full_flag = 0 
		!	smart_restart_cur_flag = 0
		!endif

		!reset TABU list
		if (TABU_flag.eq.1) then
			TABU_ban_length = 0
			TABU_allowed_length = sitenum
		endif

		call generatepdb(zero, (step+1), zero, group)
		open(5, file="backup4backbone.txt", status="replace")
			do i=1, sitenum
				write(5,704) groupdata_backup(i)%coo(1), groupdata_backup(i)%coo(2), groupdata_backup(i)%coo(3)
			enddo
			write(5,"(i1)") flag4conformer
			write(5,"(f6.3)") energy_min
			write(5,"(a6,i4)") "cycle=", SA_cycle
			write(5,"(a12, i4)") "final step=", step
		close(5)
	endif

	!*********************************
	!NOW BEGIN THE MAIN DESIGN PROCESS
	!*********************************

	!if we doing simulated annealing, then select starting temperature
	if (SA_flag.eq.1) then
		ekt_seq = SA_T_start
		ekt_backbone = SA_T_start
	endif

	do SA_cycle=1, SA_num_cycles
		!first randomize the sequence to exit out of local minimum
		if (SA_flag.eq.0.or.SA_cycle.eq.1) then !if doing simulated annealing, only want to randomize once at beginning of design
			do step=1, randomize_length
				call sequence_randomization(group, groupdata_backup)				
			enddo
			SA_best_group=group !set best simulated annealing group as randomized group (this avoids getting stuck in starting local minimum)
			!also set up the hydration properties for best peptide 
			SA_best_hyd = cur_hyd
			SA_above_hyd = above_hyd
			SA_below_hyd = below_hyd
			SA_trp_count = trp_count
		endif

		!now collect parameters for peptide and calculate/store binding energy
		call energy_parameter(group, group_para, gnum) !this loads the amber force field parameters for the pdb info into the program
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
		call bindingenergy(group, group_para, W_numex, W_inb, W_numex4, W_inb4, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
		if (energy_min.gt.1000) then
			!if energy_min has not yet been obtained, store it as the current binding energy
			energy_min = binding_energy_old
		endif 
		open(5, file="energyprofile.txt", access="append") !creates output file for energy calculations
			write(5,6) SA_cycle, 0, binding_energy_old, vdw_old, (ele_old+sgb_old), snp_old
		close(5)

		! write pdbfile
		call generatepdb(SA_cycle, zero, zero, group) 
		
		!now design peptide!
        step = 1
		do while (step.le.steps_per_cycle)
			if(step.le.steps_before_conf_change) then !!for initial first steps, don't do conformation changes since bound state is not stable (conformation changes not very helpful)
				call sequence_optimization(group, SA_best_group, groupdata_backup, step, SA_cycle, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
			else
				call ran_gen(ran2,zero)
				if(ran2.ge.seq_switch) then
					!change peptide sequence
					call sequence_optimization(group, SA_best_group, groupdata_backup, step, SA_cycle, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
				elseif (ran2.ge.backbone_switch) then
					!change peptide backbone
					call backbone_optimization(group, SA_best_group, groupdata_backup, step, SA_cycle, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old)
				else
					call rigid_body_move(group, SA_best_group, groupdata_backup, step, SA_cycle, binding_energy_old, vdw_old, ele_old, sgb_old, snp_old, group_para)
				endif
			endif

			! store data after making modification
			open(5, file="energyprofile.txt", access="append")
				write(5,6) SA_cycle, step, binding_energy_old, vdw_old, (ele_old+sgb_old), snp_old
			close(5)
			
			rmsd_flag = 1

			if (rmsd_flag.eq.1) then
				open(5, file="rmsd.txt", access="append")
					call rmsd_calculation(group, rmsd)
					rmsd_flag = 0
					write(5,703) SA_cycle, step, rmsd
				close(5)
			endif

			!ONLY UNCOMMENT BELOW SECTION FOR DEBUGGING PURPOSES
			!if (trp_count.gt.trp_limit .or. trp_count.lt.0) then
			!	open(5, file="error.txt", access="append")
			!	write(5,*) "Number of tryptophans in the peptide is outside specified range!"
			!	write(5,*) "There are currently", trp_count, "tryptophans, but the max is ", trp_limit
			!	write(5,*) "This problem is on the programmer, not the user. Terminating program."
			!	close(5)
			!	stop
			!endif

			open(5, file="energydecomp.txt", access="append")
				write(5,702) "cycle", SA_cycle, "step", step 
				do i=1,sitenum
					write(5,701) group(i)%gtype, energy_per_res(i,1), energy_per_res(i,2), energy_per_res(i,3), energy_per_res(i,4)
				enddo
			close(5)

			!periodically write the system into a pdb file
			if(mod(step,PDB_gen_freq).eq.0) then
				call generatepdb(SA_cycle, (step+1), zero, group)
				open(5, file="backup4backbone.txt", status="replace")
					do i=1, sitenum
						write(5,704) groupdata_backup(i)%coo(1), groupdata_backup(i)%coo(2), groupdata_backup(i)%coo(3)
					enddo
					write(5,"(i1)") flag4conformer
					write(5,"(f6.3)") energy_min
					write(5,"(a6,i4)") "cycle=", SA_cycle
					write(5,"(a12, i4)") "final step=", step
				close(5)
			endif

            !update dynamic cycle info
            !if (smart_restart_flag.eq.1.and.step.ge.smart_restart_wait_len) then
			!	call Smart_Restart_Update(smart_restart_cur_score, smart_restart_scores, smart_restart_index, binding_energy_old, &
			!	smart_restart_cur_len, smart_restart_cur_flag, smart_restart_best_score, smart_restart_full_flag)
			!	if (smart_restart_cur_flag.eq.1) then
			!		!current cycle not doing so hot, so do a random restart
			!		step = steps_per_cycle
			!		smart_restart_cur_len = 0
           	!		smart_restart_index = 1
            !		smart_restart_full_flag = 0 
			!		smart_restart_cur_flag = 0
			!		if (net_step.lt.nstep_terminal.and.cycle.ge.num_restarts) then
			!			!keep design going at least until nsteps are complete if using smart subcycle
			!			cycle = cycle - 1
			!		endif
			!	endif
            !endif

			!update residue selection probability matrix, if applicable
			!if (decomp_mutate_flag.eq.1) then
			!	call calc_mutate_ps(energy_per_res, res_select_p_vector)
			!endif
			net_step = net_step + 1
			step = step + 1
		enddo

		!do simulated annealing bookkeeping
		if (SA_flag.eq.1) then
			ekt_seq = ekt_seq * SA_fraction
			ekt_backbone = ekt_backbone * SA_fraction
			group = SA_best_group
			cur_hyd = SA_best_hyd 
			above_hyd = SA_above_hyd
			below_hyd = SA_below_hyd
			trp_count = SA_trp_count
			open(5, file="output.txt", access="append")
			write(5,*) "Starting new cycle, reference kBT value is ", ekt_seq
			close(5)
		endif

		!ONLY UNCOMMENT BELOW SECTION FOR DEBUGGING PURPOSES
		!check the hydration properties are still okay
		!do i=1,6
		!	if (cur_hyd(i).gt.max_hyd(i) .or. cur_hyd(i) .lt. min_hyd(i) ) then
		!			open(5, file="error.txt", access="append")
		!		write(5,*) "Bad hydration property found!"
		!		write(5,*) "For index ", i, "value is ", cur_hyd(i)
		!		write(5,*) "For index ", i, " the minimum is ", min_hyd(i), " and the maximum is ", max_hyd(i)
		!		close(5)
		!		stop
		!	endif
		!enddo
	enddo
6	format(i5, i5, 4f10.4)
701 format(a5, 4f10.3)
702 format(a6, i4, a7, i5)
703 format(i3, i5, f7.3)
704 format(3f6.3)

	!write out probability data
	open(102, file='output.txt', access='append')
	write(102,*); write(102,*) "********************"
	write(102, *) "Probalities of Accepting Peptide Changes"
	if (SEQ_SCMF_attempts.eq.0) SEQ_SCMF_attempts = 1
	write(102,705) "SCMF Mutations", SEQ_SCMF_success*1./SEQ_SCMF_attempts
	if (SEQ_POINT_attempts.eq.0) SEQ_POINT_attempts = 1
	write(102,705) "Point Mutations", SEQ_POINT_success*1./SEQ_POINT_attempts
	if (CONF_N2C_attempts.eq.0) CONF_N2C_attempts = 1
	write(102,705) "N and C termini CONROT moves", CONF_N2C_success*1./CONF_N2C_attempts
	if (CONF_N_attempts.eq.0) CONF_N_attempts = 1
	write(102,705) "N terminus CONROT moves", CONF_N_success*1./CONF_N_attempts
	if (CONF_C_attempts.eq.0) CONF_C_attempts = 1
	write(102,705) "C terminus CONROT moves", CONF_C_success*1./CONF_C_attempts
	if (CONF_Internal_attempts.eq.0) CONF_Internal_attempts = 1
	write(102,705) "Internal CONROT moves", CONF_Internal_success*1./CONF_Internal_attempts
	write(102, *); write(102,*) "Per Residue probabilities of successful sequence changes"
	do i=1, sitenum
		if (SEQ_PerRes_attempts(i).eq.0) SEQ_PerRes_attempts(i) = 1
		write(102,706) "Residue", i, SEQ_PerRes_success(i)*1./SEQ_PerRes_attempts(i)
	enddo
	close(102)
	705 format(a30, f8.3)
	706 format(a8, i3, f8.3)

	!calculate the non-polar solvation energy for the top 50 peptides, if desired.
	!if (ESURF_flag.eq.1) then
	!	open(101, file="ESURF_Results.txt", access='append')
	!	write(101,*) "ORIGINAL    ESURF      NET"
	!	do i=1,ESURF_count
	!		if (top_pep_scores(i) > 100) exit !ESURF array not completely filled up and rest is also empty, so end program
	!		
	!		!read the pdb file
	!		fname = "pdbfiles/"//ESURF_pep_descriptors(i)
	!		call readpdb_for_ESURF(group, gnum, fname)
	!			
	!		!load the force field parameters for system (note that atom links don't need to be built here, not used in calculation)
	!		call energy_parameter(group, group_para, gnum) !this loads the amber force field parameters for the pdb info into the program
	!
	!		!non-polar solvation energy for whole system
	!		call arvo(one, gnum, group, group_para, volume, surfarea)
	!		comp_ESURF=surftens*surfarea+offsetvalue
	!
	!		!non-polar solvation energy for receptor
	!		call arvo(sitenum+1, gnum, group, group_para, volume, surfarea)
	!		rec_ESURF=surftens*surfarea+offsetvalue
	!
	!		!non-polar solvation energy for peptide
	!		call arvo(one, sitenum, group, group_para, volume, surfarea)
	!		lig_ESURF=surftens*surfarea+offsetvalue
	!
	!		!delta non-polar solvation energy
	!		net_esurf = comp_ESURF - rec_ESURF - lig_ESURF
	!		write(101,*) (group(step)%gtype, step=1, sitenum) 
	!		write(101, *) ESURF_pep_descriptors(i)
	!		write(101, 707) top_pep_scores(i), net_ESURF, top_pep_scores(i)+net_ESURF
	!		write(101,*)
	!	enddo
	!	close(101)
	!endif
707	format(3f9.3)

    !if (smart_restart_flag.eq.1) then
    !    deallocate(smart_restart_scores)
    !endif
end