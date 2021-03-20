# Recreates all the visualizations of the model output for the paper

WORKING_DIR=$(PWD)
DESIRED_DIR="capability-accumulation"

if  [ ${WORKING_DIR[-23,-1]} != $DESIRED_DIR ]; then
    cd ".."
    WORKING_DIR=$(PWD)
    if  [ ${WORKING_DIR[-23,-1]}==$DESIRED_DIR ]; then
        echo "Working directory was changed to:" $WORKING_DIR
    else 
        echo "Still wrong working directory:" $WORKING_DIR
        exit 1 
    fi
else
    echo "Working directory correct:" $WORKING_DIR
fi
  
echo "Start visualization"
python python/convert_tmax.py
Rscript R/vis_topology_dynamics.R
Rscript R/vis_topology_finaldist.R
Rscript R/vis_allocation_finaldist.R
Rscript R/vis_sensitivity_finaldist.R
Rscript R/vis_topdist_finaldist.R
Rscript R/vis_vistop_finaldist.R
echo "Completed visualization"