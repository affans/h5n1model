# ABM for H5N1 Transmission Dynamics

## Environment
- Julia 1.11.5 
- Use `Pkg` to `resolve` and `instantiate` the environment in `Project.toml`. 

## Running the model
- Start Julia in the supplied environment.
- Include the model by `include("model.jl")`.
- The function `time_loop()` runs the model (i.e., a single realization) with the supplied parameters. Parameters (with default values) are 
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
> The results in the manuscript are based on a `n` independent runs of `time_loop()` for each scenario (see script.jl file). This was done in a parallel fashion on a compute cluster (see `script.jl`). One may also just use a for loop to run n replicates, but this might take some time. Each independent simulation is only a found seconds though. 

> [!NOTE] 
> We provide output of all scenarios and simulations. The data can be found in a related repository found here: 
[Data Repository](http://github.com/affans/h5n1model_data).


- The function `calibrate()` is used to calibrate and find the number of secondary infections. The function returns the basic reproduction number. The parameters include number of simulations, beta value, and isolation parameters. Calibration involved finding a beta value such that the number of secondary infections equals some reproduction number. The

> [!NOTE]  
> In our manuscript, we provide results We have already calibrated the model to a reproduction The beta value `0.017` corresponds to a reproduction number of 0.2.  

## Individual Modules

To test the model's individual modules, one can `include("model.jl")` and start the initialization process. The following initialization functions are needed before the model is run: 
  - `init_state()` initializes the model variables. 
  - `init_agents()` initializes agents and sets up demographics.
  - `init_households2()` initializes household composition and associates agents to each household
  - `init_farming()` initializes farms and associates agents to farms as workers.   

> [!NOTE]  
> A helper function `init_model()` runs these functions in the correct order

> [!NOTE]  
> To enable debugging, set up a debug logger `using Logging; debuglogger = ConsoleLogger(stderr, Logging.Debug); global_logger(debuglogger)`

After initialization, an infected person can be introduced by `insert_infection()`. The natural history can be debugged using `natural_history()` and `activate_swaps()` functions (they will require the agent as an input parameter). The function `transmission_with_contacts()` (with agent passed in as the parameter) determines the number of contacts, with who they will contact, and whether disease transmission will occur. The function `virus_mutation()` calculates the probability of virus mutation for each agent specifically. 

## The `script.jl` file

The `script.jl` file is an internal file not required for the model. It was used to calibrate the model, set up the scenarios described in the paper,  activates a SLURM HPC cluster, and submits simulations in a parallel fashion. Results are collected, processed, and saved to disk for further analysis. The file is provided as a template for parallel simulations and processing of simulation data. 

