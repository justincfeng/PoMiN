############################################################
# test.sh - test known cases for correct results
############################################################
#!/bin/sh
make -s quad conv

############################################################
# CONVERGENCE
############################################################
read -r -d '' input << 'EOF'
orbit-0.00 roughly modelled after Mercury and the Sun,0,1507000000000000,15070000000000,0,0,1,1,1,0,0,0,0,-2.6716223E-12,0,0.000000166,3860700000,0,0,0,2.6716223E-12,0
orbit-0.50 newtonian orbit with an eccentricity of 0.50,0,4.26308E+015,42630800000000,0,0,1,1,1,0,0,0,0,-3.2720557E-12,0,0.000000166,3860700000,0,0,0,3.2720557E-12,0
orbit-0.90 newtonian orbit with an eccentricity of 0.90,0,4.76627E+016,47662700000000,0,0,1,1,1,0,0,0,0,-3.6825772E-12,0,0.000000166,3860700000,0,0,0,3.6825772E-12,0
orbit-0.95 newtonian orbit with an eccentricity of 0.95,0,1.34811E+017,135000000000000,0,0,1,1,1,0,0,0,0,-3.7307175E-12,0,0.000000166,3860700000,0,0,0,3.7307175E-12,0
scatter-massless-massive massless particle deflecting off stationary massive particle,0,100000,1000,0,0.1,1,1,1,0,0,0,-2E-12,0,0,0,-50000,500,0,2E-12,0,0
scatter-massless-massive massless particle deflecting off stationary massive particle,0,250,2.5,0,0,1,1,1,-0.0000001,2.82048732949973E-12,0,-2.00796147822064E-12,7.75749609071621E-15,0,0,-19.5519457042,498.093058565,0,2.00796147822064E-12,-7.75749609071621E-15,0
scatter-massive-massive massive particle deflecting off stationary massive particle,0,100000,1000,0,0,1,1,60,0,0,0,-0.000002,0,0,0.0000002,-50000,1000,0,0.000002,0,0
scatter-massless-massless deflection of two photons,0,100000,1000,0,0.1,1,1,0,50000,50,0,-0.001,0,0,0,-50000,-50,0,0.001,0,0
scatter-massless-massless deflection of two photons,0,25,0.25,0,0,1,1,0,1.0677612488,49.99848548,0,-0.001000031,-3.10637499649484E-08,0,0,-1.06776125,-49.9984854741,0,0.001000031,3.10637499649484E-08,0
EOF

while read -r line; do
    # Divide the timestep value by $i and run that as input.
    for (( i=1; i<=512; i=$i*2 )); do
        timestep=$(echo $line | awk -F, '{print $4}')
        timestep=$(echo $timestep $i | awk '{printf("%.30e", $1/$2)}')
        
        # Append the last line.
        last_lines+="$(printf "$line" | awk -F, -v OFS=, '{$4=ts; print }' ts=$timestep \
            | ./pomin $@ | tail -n1)\n"
    done

    # Print header for readability.
    printf "\n${line%%,*}\n"
    # Calculate convergence factors.
    printf "$last_lines" > temp
    ./validation/conv/convergence.exe temp

done <<< "$input"

############################################################
# CLEANUP
############################################################
# Delete executables.
rm ./validation/conv/convergence.exe ./pomin

exit 0
############################################################
# END test.sh
############################################################
