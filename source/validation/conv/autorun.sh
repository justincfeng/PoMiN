############################################################
# autorun.sh - consecutively run pomin(1) inputs 
############################################################
#!/bin/bash

# Compile code and convergence.c script.
make quad conv

# Define Initial Filenames
FILECF=1
FILECV=2

# Read from stdin.
while read N a line; do

    # Total number of particles to read
    #    N=5 
        # $(echo $line | awk -F, '{print $2+0}')

    # Print header 
    printf "iteration,time" > $FILECF
    printf "iteration,time" > $FILECV
    for (( i=1; i<=$N; i++ )); do
        printf ',qx_%s,qy_%s,qz_%s,px_%s,py_%s,pz_%s'\
            "$i" "$i" "$i" "$i" "$i" "$i" >> $FILECF
    if [ $i = $a ]; then
        printf ',qx_%s,qy_%s,qz_%s,px_%s,py_%s,pz_%s'\
            "$i" "$i" "$i" "$i" "$i" "$i" >> $FILECV
    fi

    done
    printf "\n" >> $FILECF
    printf "\n" >> $FILECV

    # Half the timestep value and run that as input.
    for (( i=1; i<=512; i=$i*2 )); do
        timestep=$(echo $line | awk -F, '{print $4}')
        timestep=$(echo $timestep $i | awk '{printf("%.30e", $1/$2)}')
        
        echo $line | awk -F, -v OFS=, '{$4=ts; print }' ts=$timestep \
         | ./pomin | tail -n1 >> $FILECF
	echo $i

    fileline=`tail -n 1 $FILECF`

    k=$((2+6*(a-1)))

    cvline=$(echo $fileline | awk -F, -v x=$k '{for(j=x+1;j<=x+6;j++){printf ",%s", $j}; printf "\n"}')

    echo $i,$timestep$cvline >> $FILECV

    done

    # Print header for readability.
    printf "\nQ factors:\n"
    # Calculate convergence factors.
    ./validation/conv/convergence.exe $FILECV

rm ./1
rm ./2

done

############################################################
# CLEANUP
############################################################
# Delete executables.
rm ./validation/conv/convergence.exe #./pomin

exit 0
############################################################
# END autorun.sh
############################################################

