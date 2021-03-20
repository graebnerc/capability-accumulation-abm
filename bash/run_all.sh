
# Runs the following commands in parallel:
# zsh bash/run_allocation_comp.sh
# zsh bash/run_comp_dist.sh
# zsh bash/run_deltas.sh
# zsh bash/run_saturation.sh
# zsh bash/run_topologies.sh
# zsh bash/run_topdist.sh
# zsh bash/run_vistop.sh
# zsh bash/run_range_of_vision.sh

echo "Start parallel session"
parallel --verbose --ungroup zsh ::: bash/run_allocation_comp.sh bash/run_comp_dist.sh bash/run_deltas.sh bash/run_range_of_vision.sh bash/run_topologies.sh bash/run_saturation.sh bash/run_topdist.sh bash/run_vistop.sh 
echo "Completed parallel session"

python python/convert_tmax.py

Rscript R/visualize.R