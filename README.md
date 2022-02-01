## Fish CPR task
Utility code for the research article
"Mine or ours? Neural basis of the exploitation of common-pool resources" by
Mario Martinez-Saito, Sandra Andraszewicz, Vasily Klucharev, Jörg Rieskamp.

MRI data processing (Matlab with SPM) and data analysis (Matlab) files. 

### Index
- Main file
  * sazan_main.m: behavioral model fitting, analysis, and plotting 
- File reading routines
  * sazan_readlog.m
- MRI preprocessing with SPM 12 and analysis routines
  * mybatch_sazan.m
- Preprocessed data
  * sazan_fitted.mat
- Data analysis routines
  * sazan_nll.m: negative log-likelihood function
  * sazan_optim.m: model optimization
  * sazan_sim.m: simulations for posterior predictive analysis
- License
  * LICENSE: MIT license
