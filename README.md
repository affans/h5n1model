# ABM for H5N1 Transmission Dynamics

## Requires
- Julia 1.11.5 
- Use `Pkg` to `resolve` and `instantiate` the environment in `Project.toml`. 

## Running the model
- Start Julia and enable the environment.
- Include the model by `include("model.jl")`.
- The function `time_loop(...)` runs the model (i.e., a single realization) with the supplied parameters. Parameters (with default values) are 
  ```
    simid=1,     # simulation ID
    beta=0.0,    # probability of transmission 
    init_inf=0,  # initial infected individuals
    iso_day=-1,  # will isolate after this many days
    iso_prop=0.0 # proportion of islation
    vac_scn1=NONE,  # vax scenario: NONE, FARM, FARM/HH
    vac_scn2=A0, # vax type: against infection, symptomatic, or both
    vac_cov=0.0, # vaccine coverage
    vac_time=42  # when vaccine will start
  ```
  > [!NOTE]  
  > The function will automatically initialize the model (so no need to run `init_` functions - see below)

  > [!NOTE]  
  > The results in the manuscript are based on a `n = 1000` independent runs of `time_loop()` for each scenario (about 500 scenarios). This was done in a parallel fashion on a compute cluster (see `script.jl`). One may also just use a for loop to run 1000 replicates. 

- The function `calibrate(...)` is used to calibrate and find the number of secondary infections. The parameters include number of simulations, beta value, and isolation parameters. Calibration involved finding a beta value such that the number of secondary infections equal 1.2, 1.5, and 1.8 (reproduction numbers). Isolation was turned on after calibration to see the reduction in secondary infections due to isolation.
  > [!NOTE]  
  > The results in the manuscript are based on a `n = 1000` independent runs of `calibrate` for each scenario. This was done in a parallel fashion on a compute cluster (see `script.jl`). One may also just use a for loop to run 1000 replicates. 

## Components of the model

To test the model's components, one can `include("model.jl")` and start the initialization process. The following initialization functions are needed before the model is run: 
  - `init_state()` initializes the model variables. 
  - `init_agents()` initializes agents and sets up demographics.
  - `init_households2()` initializes household composition and associates agents to each household
  - `init_farming()` initializes farms and associates agents to farms as workers.   

> [!NOTE]  
> A helper function `init_model()` runs these functions in the correct order

> [!NOTE]  
> To enable debugging, set up a debug logger `using Logging; debuglogger = ConsoleLogger(stderr, Logging.Debug); global_logger(debuglogger)`

After initialization, an infected person can be introduced by `insert_infection()`. The natural history can be debugged using `natural_history()` and `activate_swaps()` functions (they will require the agent as an input parameter). The function `transmission_with_contacts()` (with agent passed in as the parameter) determines the number of contacts, with who they will contact, and whether disease transmission will occur. 

## The `script.jl` file
This file runs the model for specified scenarios (as presented in the table). The model activates a SLURM HPC cluster and submits simulations in a parallel fashion. Results are collected, processed, and saved to disk. It's unlikely that `script.jl` can be used directly for reproducing the results.