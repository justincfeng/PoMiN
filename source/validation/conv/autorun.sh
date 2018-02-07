############################################################
# autorun.sh - consecutively run pomin(1) inputs 
############################################################
#!/bin/bash

# Compile code and convergence.c script.
make quad conv

# Define an iterator over particle number.
a=1

# Read from stdin.
while read line; do

    # Read number of particles.
    N=2 # $(echo $line | awk -F, '{print $2+0}')

    # Print header 
    printf "iteration,time" > $a
    for (( i=1; i<=$N; i++ )); do
        printf ',qx_%s,qy_%s,qz_%s,px_%s,py_%s,pz_%s'\
            "$i" "$i" "$i" "$i" "$i" "$i" >> $a
    done
    printf "\n" >> $a

    # Half the timestep value and run that as input.
    for (( i=1; i<=512; i=$i*2 )); do
        timestep=$(echo $line | awk -F, '{print $4}')
        timestep=$(echo $timestep $i | awk '{printf("%.30e", $1/$2)}')
        
        echo $line | awk -F, -v OFS=, '{$4=ts; print }' ts=$timestep \
         | ./pomin | tail -n1 >> $a
	echo $i
    done

    # Print header for readability.
    printf "\nQ factors:\n"
    # Calculate convergence factors.
    ./validation/conv/convergence.exe $a

    ((a = a + 1))
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

