# met-uncertainty-paper

## Instructions

Download `daymet_v4_stnxval_tmax_na_2022.nc` from [ORNL DAAC](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2132) and copy into `data` folder.

1. Run

```
julia loo_validation.jl -t [number of threads]
```

to do the main leave-one-out cross-validation and write the results to
`output`. This takes ~20 min using 4 threads on my laptop.

2. Run 

```
julia agg_loo_res.jl -t [number of threads]
```

to aggregate the results from the previous step and compute the UQ accuracy
statistics using nested symmetric prediction intervals.

3. `plot_loo_valid_res.jl` creates figures 1 and 2\*.

4. `full_fig.jl` creates figure 3\*.

5. `krig_vs_sim.jl` compares simulations from kriging versus conditional
simulations for two nearby stations. It also creates figures 4 and 5\*.

\* The scripts as written here use `GLMakie` for plotting, which opens a local
interactive viewing window. To save the figures as PDFs, uncomment
`using CairoMakie` and the the `save(...)` lines in each file.
