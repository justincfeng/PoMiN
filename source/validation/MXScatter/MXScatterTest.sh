############################################################
# MXScatterTest.sh - Momentum Exchange Scattering Test
# 
# Script assumes initial data in CM frame and G=c=1
############################################################

# Compile pomin
# make 

# Read from stdin.
while read line; do

	# Pick out 1/2 the impact parameter (y-value of particle 2)
	hb=$(echo $line | awk -F, '{print $19}')

	# Pick out initial y-momentum for particle 1
	pstrt=$(echo $line | awk -F, '{print $15}')

	# Pick out mass of particle 1
	m1=$(echo $line | awk -F, '{print $10}')

	# Pick out x-momentum of particle 1
	p=$(echo $line | awk -F, '{print $14}')

	# Pick out mass of particle 2
	m2=$(echo $line | awk -F, '{print $17}')

	# Run pomin and get last line of output
    EndL=$(echo $line | ./pomin | tail -n1)

	# Final y-momentum of particle 1
	pfin=$(echo $EndL | awk -F, '{print $7}')

	# Compute impact parameter
	b=$(echo $hb | awk -v PREC="quad" '{bi=2*$1; printf("%0.40e", bi)}')

	# Compute energy of particle 1
	E1=$(echo $m1 $p | awk -v PREC="quad" '{en = sqrt(($1)^2+($2)^2); printf("%0.40e", en)}')

	# Compute energy of particle 2
	E2=$(echo $m2 $p | awk -v PREC="quad" '{en = sqrt(($1)^2+($2)^2); printf("%0.40e", en)}')

	# Compute theoretical change in momentum
    pDT=$(echo $hb $p $E1 $E2 | awk -v PREC="quad" '{pd = ($3**2*$4**2*(1 + (1/($3**2) + 1/($4**2) + 4/($3*$4))*$2**2 + $2**4/($3**2*$4**2)))/($1*($3 + $4)*$2); printf("%0.40e", pd)}')

	# Compute difference between initial and final y-momentum
	pDIFF=$(echo $pfin $pstrt | awk -v PREC="quad" '{diff = $1-$2; printf("%0.40e", diff)}')

	# Compute % difference between analytical and numerical result
	Res=$(echo $pDT $pDIFF | awk -v PREC="quad" '{R = 100*(sqrt(($1-$2)^2)/$1); printf("%0.40e", R)}')

	# Send results to stdout	
	echo $b,$pDT,$pDIFF,$Res

done

############################################################
# CLEANUP
############################################################
# Delete executables.
# rm ./pomin

exit 0
############################################################
# END MXScatterTest.sh
############################################################

