package require pbctools
 #--------------------------------------------------------------------------------------------#
#VMD Stuff

set step 1

#if using the dehydrated trajectories
set in_dmsfile [pwd]/[glob *dehydrated.dms]; #pwd gives the path and .dms gives the dms file
set dcd_file [pwd]/[glob *dehydrated.dcd]
mol new $in_dmsfile
mol addfile $dcd_file type dcd first 0 last -1 step $step waitfor all

#if using original trajectory
#set dtr_clickme [pwd]/[glob *_trj/clickme.dtr]; #pwd gives the path and .dms gives the dms file
#set in_cmsfile [pwd]/[glob *-in.cms]
#mol new $in_cmsfile
#mol addfile $dtr_clickme type dtr first 1 last -1 step $step waitfor all

#--------------------------------------------------------------------------------------------#
#Global Variables
set nf [molinfo top get numframes]
set intv 0.0096 ; #in nanoseconds
##--------------------------------------------------------------------------------------------#
set fh1 [open "cv_results.dat" w]
puts $fh1 "# this data file was made with $argv"
puts $fh1 "# Chains E,F,G,H are on residues 50-69 of chains A through D (present only in the C8 cage models!)"
puts $fh1 "# Chains I,J,K,L are on residues 10-29 of chains A through D (present in both C4 and C8 cage models)"
puts $fh1 "# this script generates CV's for the shotgun simulation restarts"
puts $fh1 "# \[frame, distance along z cordinate of protein (0 is center), distance protein is away from the center of the top of the cage in the plane]"
#puts $fh1 "#--------------------------------------------------------------------------------------------"

#--------------------------------------------------------------------------------------------#

#set chain_list [lsort -unique [[atomselect top "resname D1D and noh" frame 0] get {chain}]]
#puts $fh1 [format "!CHAIN %s" $chain_list]

set prot [atomselect top "chain P and noh" frame 0]
set cage_top [atomselect top "(chain A to D) and (resid 10 to 29) and noh" frame 0]
set ca [atomselect top "(chain A) and (resid 17 to 22) and noh" frame 0]
set cb [atomselect top "(chain B) and (resid 17 to 22) and noh" frame 0]
set cc [atomselect top "(chain C) and (resid 17 to 22) and noh" frame 0]
set cd [atomselect top "(chain D) and (resid 17 to 22) and noh" frame 0]
set full_cage [atomselect top "(chain A to D) and noh" frame 0]

for {set i 0} {$i<$nf} {incr i} {

	set cvs {}
	lappend cvs $i
	animate goto $i
	$prot frame $i
	$cage_top frame $i
	$full_cage frame $i
	$ca frame $i
	$cb frame $i
	$cc frame $i
	$cd frame $i
	pbc unwrap -sel "all and not element Mg and not element Na and not water" -now

	#CV 1: projection of cage to protein vector onto cage to cage top vector
	set hsa_xyz [measure center $prot]
	set cagetop_xyz [measure center $cage_top]
	set cage_xyz [measure center $full_cage]
	set v_cage_cagetop [vecsub $cagetop_xyz $cage_xyz]; #gives vector from com of cage to com of top of cage (chains I,J,K,L)
	set norm [vecnorm $v_cage_cagetop]; #normalizes vector
	set v_cage_prot [vecsub $hsa_xyz $cage_xyz]; #vector of cage center to protein
	set cv1 [vecdot $v_cage_prot $norm]


	#CV2: projection of protein to the top of chain A (face of cage), onto center of the top of the cage vector to the side of the cage vector (makes sure it doesnt stray over the edge.

	#xyz coords for the center point of all four edges (chain a-d)
	set ma [measure center $ca]
	set mb [measure center $cb]
	set mc [measure center $cc]
	set md [measure center $cd]

	#vector from hsa to the side of the cage
	set v_hsa [vecsub $hsa_xyz $cagetop_xyz]


	#create two vectors the disect the top of the cage
	set vec1 [vecsub $ma $mc]
	set vec2 [vecsub $mb $md]
	#normalize them
	set nvec1 [vecnorm $vec1]
	set nvec2 [vecnorm $vec2]
	#cross product gives an orthagonal vector
	set ncvec [vecnorm [veccross $nvec1 $nvec2]]
	#second cross product creates a truely orthagonal vector which defines the plane (better than the original two vectors)
	set cvec2 [veccross $ncvec $nvec1]

	#puts [veclength $cvec2]
	#puts [vecdot $cvec2 $nvec1]

	set proj_1 [vecdot $v_hsa $nvec1]
	set proj_2 [vecdot $v_hsa $cvec2]
	#set proj_3 [vecdot $v_hsa $ncvec] #dont need because we have the vector i already defined
	set test {}
	lappend test $proj_1
	lappend test $proj_2
	set dist_from_center [veclength $test] ;#this is the distance of the protein from the center of the top of the cage as defined by the plane given in nvec1 and cvec2
	puts $fh1 [format "%i,%f,%f" $i $cv1 $dist_from_center]
	flush $fh1

}


#--------------------------------------------------------------------------------------------#
exit


"""
	#vector from the top face center of the cage to the side of the cage
	set v_cent_ma [vecsub $ma $cagetop_xyz]
	set v_cent_mb [vecsub $mb $cagetop_xyz]
	set v_cent_mc [vecsub $mc $cagetop_xyz]
	set v_cent_md [vecsub $md $cagetop_xyz]

	#normaled cagetop to side vector in four directions
	set norm_sa [vecnorm $v_cent_ma]
	set norm_sb [vecnorm $v_cent_mb]
	set norm_sc [vecnorm $v_cent_mc]
	set norm_sd [vecnorm $v_cent_md]

	#dot product for magnitude
	set cv2 [vecdot $v_hsa_ma $norm_sa]
	set cv3 [vecdot $v_hsa_mb $norm_sb]
	set cv4 [vecdot $v_hsa_mc $norm_sc]
	set cv5 [vecdot $v_hsa_md $norm_sd]

	#puts $fh1 [format "%i , %.4f , %.4f, %.4f" $i $d1 $d2 $magnitude]
	#lappend cvs $cv1
	#lappend cvs $cv2
	#lappend cvs $cv3
	#lappend cvs $cv4
	#lappend cvs $cv5
	#puts $fh1 $cvs
	#puts $cvs
"""
