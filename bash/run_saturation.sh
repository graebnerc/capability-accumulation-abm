# This file runs the simulation and then calls the visualization script

# Set the parameters here--------------
IDENTIFIER="saturation"
RUNS="50"
YDIM="all"
COLVAR="nominal_demand"

# The rest of the script should not be altered-------------

echo "----------------------------------------------------"
echo "Run MCS"
echo "----------------------------------------------------"

FULLID="python/parameters/"$IDENTIFIER
AGGFILE="output/"$IDENTIFIER"/"$IDENTIFIER"_agg.feather"

simulate(){
python python/MCS.py $FULLID $RUNS
}
simulate

echo "----------------------------------------------------"
echo "Finished!"
echo "----------------------------------------------------"
