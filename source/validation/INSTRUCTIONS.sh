################################################################################
# CONVERGENCE TESTS
################################################################################

# If you run this file as a script, run it in this directory using the command: 
# bash INSTRUCTIONS.sh

cd ..

# To run the convergence tests, run the following from the source directory:

a1=`grep orbit-0.50 input/input.csv`
echo 2 2 "$a1" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_Orbit-050.txt

a2=`grep orbit-0.90 input/input.csv`
echo 2 2 "$a2" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_Orbit-090.txt

a3=`grep orbit-0.95 input/input.csv`
echo 2 2 "$a3" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_Orbit-095.txt

a4=`grep SCMS1- input/input.csv`
echo 2 2 "$a4" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_SCMS1.txt

a5=`grep SCML1- input/input.csv`
echo 2 2 "$a5" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_SCML1.txt

a6=`grep SCMX1- input/input.csv`
echo 2 2 "$a6" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_SCMX1.txt

a7=`grep CFiveParticles input/input.csv`
echo 5 2 "$a7" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_Five.txt

a8=`grep TPFiveParticles input/input.csv` 
echo 5 4 "$a8" | bash validation/conv/autorun.sh > validation/conv/results/raw/conv_TPFive.txt

# The outputs may be found in the directory: validation/conv/results/raw/

# Each variable aN (Ex: a7) corresponds to a line of initial data found in the
# directory input/input.csv

# autorun.sh takes inputs of the form: n1 n2 "$aN" 
#  n1 is total number of particles
#  n2 is the particle on which to apply the convergence test 
#   (convergence.c is set up to apply the convergence test to the quantity p^2)

# NOTE: For each variable aN (Ex: a7), make sure grep returns one line only; you 
#       will get errors if grep returns multiple lines. To test this, run the 
#		command indicated in each variable. For example, to test a7, run the
#		command:
#		grep CFiveParticles input/input.csv



################################################################################
# MOMENTUM EXCHANGE SCATTERING TESTS
################################################################################

# To run the scattering tests, first run make in the source directory:

make

# Next, run the following command lines from the source directory:

grep CMScatter_b input/MXIC/inputscatter_b_massive.csv | bash validation/MXScatter/MXScatterTest.sh > validation/MXScatter/results/raw/scatter_b_massive.csv

grep CMScatter_b input/MXIC/inputscatter_b_massless.csv | bash validation/MXScatter/MXScatterTest.sh > validation/MXScatter/results/raw/scatter_b_massless.csv

grep CMScatter_b input/MXIC/inputscatter_b_mixed.csv | bash validation/MXScatter/MXScatterTest.sh > validation/MXScatter/results/raw/scatter_b_mixed.csv

# To insert a header into these files, run the following:

echo 'Impact Parameter,Analytical dp,Numerical dp,Percent difference' | cat - validation/MXScatter/results/raw/scatter_b_massive.csv > temp && mv temp validation/MXScatter/results/raw/scatter_b_massive.csv

echo 'Impact Parameter,Analytical dp,Numerical dp,Percent difference' | cat - validation/MXScatter/results/raw/scatter_b_massless.csv > temp && mv temp validation/MXScatter/results/raw/scatter_b_massless.csv

echo 'Impact Parameter,Analytical dp,Numerical dp,Percent difference' | cat - validation/MXScatter/results/raw/scatter_b_mixed.csv > temp && mv temp validation/MXScatter/results/raw/scatter_b_mixed.csv

# The outputs may be found in the directory: validation/MXScatter/results/raw/



############################################################
# CLEANUP
############################################################
# Delete executables.
rm ./pomin

exit 0
############################################################
# END INSTRUCTIONS
############################################################

rm ./pomin

