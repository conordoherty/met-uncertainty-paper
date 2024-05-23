# met-uncertainty-paper

## Data preparation

Download `daymet_v4_stnxval_tmax_na_2022.nc` from [ORNL DAAC](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2132) and copy into `data` folder.

## Running the code 

### Setup

In the Julia REPL, run the following commands to activate the environment
and install the dependencies.

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Predictive distribution validation for CA in 2022

There are two scripts that can be run either from the REPL using `include("[...].jl")` or from the
command line using `julia --project=. [...].jl`. Both scripts can be run multithreaded. Pass
`-t [number of threads]` when starting the REPL or running the script from the command line.

1. `loo_validation.jl` performs leave-one-out (LOO) cross-validation for each stations with
data for each day of 2022. Results for each day are saved in `./output`. This runs in ~20 min
using 4 threads on my laptop.

2. `agg_loo_res.jl` processes the output of the previous script to compute the actual coverage
of theoretical predictive distributions. Results are saved in `./output`.

### Figures

1. `plot_loo_valid_res.jl` creates figures 4 and 5. It requires having run the LOO validation
code described above.

2. `full_fig.jl` generates uncertainty maps for four different days. It creates figure 6.

3. `krig_vs_sim.jl` compares simulations from kriging versus conditional
simulations for two nearby stations. It creates figures 7 and 8.

4. `station_fig.jl` creates figure 1, which shows the locations of weather stations.

5. `valid_process_fig.jl` creates figure 2, which shows how the actual coverage of
theoretical prediction intervals is calculated.
