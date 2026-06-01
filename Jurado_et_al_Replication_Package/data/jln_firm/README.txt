The steps for constructing firm level uncertainty is similar to the
construction of macro uncertainty.
1. Collect firm level data for profits in firmpanel.xlsx (from Compustat
as  in Bloom (2009)
2. Transform data as necessary
3. Generate_ferrors.m
      Estimate factors from  firms' profits
      Compute forecast errors for each of the factors using AR(4) model
      Compute forecast errors for firm profit series using factor-augmented regression
      output: vyt.txt and vft.txt
4. Estimate stochastic volatility in the forecast errors using 'stochvol' (R package):
     generate_svfdraws.R -> svfmeans.txt
     generate_svydraws.R -> svymeans.txt
5. Use svfmeans.txt and svymenas.txt as input to individual uncertainty estimates
    generate_ut.m -> ut.mat
6. Aggregate ut into  uncertainty common to firms
    generate_firmu.m -> firmu.mat (also includes individual  uncertainty 
           estimates for each firm, and for h=1..12)

plot_firmu.m replicates Figure 8