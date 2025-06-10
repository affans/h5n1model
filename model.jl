using StatsBase, Random, Distributions, DelimitedFiles
using StaticArrays

# make sure the order doesn't change
@enum INFTYPE UNDEF=0 SUS=1 INF=2 REC=3 

# define an agent and all agent properties
Base.@kwdef mutable struct Human
    idx::Int64 = 0
    age::Int64 = 0 # age in weeks
    agegroup::Int8 = 0 # store the age group for easier access to other properties
    hid::Int64 = 0 # household id
    fid::Int64 = 0 # farm id 
    inf::INFTYPE = SUS # infection status
    swap::INFTYPE = UNDEF # swap state
    tis::Int16 = 0 # time in state
    st::Int16 = typemax(Int16) # max time in state 
    vac::Int16 = 0 # dose count
    eff::Float64 = 0.0 # vaccine efficacy
    exp::Int64 = 0 # waning time of vaccine
end
Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## system parameters
Base.@kwdef mutable struct ModelParameters
    flag_init = false # flag to indicate if the model is initialized
    flag_agents = false # flag to indicate if agents are initialized
    flag_households = false # flag to indicate if households are initialized 
    flag_farming = false # flag to indicate if farming is initialized    
end

include("helpers.jl") # include the helper functions
include("debug.jl")

# constant variables
const POPSIZE = 153565
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters

# age distribution of 0 - 100 years old
const _POPDISTR = @SVector [1530, 1545, 1565, 1670, 1700, 1820, 1755, 1790, 1820, 1935, 1875, 1925, 
                            1925, 1910, 1890, 1800, 1895, 1900, 1820, 1855, 1995, 2135, 2040, 2060, 
                            2030, 2120, 1985, 2020, 2125, 2140, 2120, 2230, 2240, 2065, 2105, 2190, 
                            2160, 2195, 2225, 2210, 2190, 2025, 2005, 1965, 1845, 1870, 1920, 1840, 
                            1835, 1855, 1880, 1895, 1870, 1915, 1645, 1820, 1885, 1970, 1875, 1755, 
                            1850, 1945, 1845, 1850, 1790, 1720, 1810, 1735, 1595, 1520, 1460, 1515, 
                            1430, 1500, 1410, 1135, 1060, 1050, 925, 860, 820, 735, 730, 755, 570, 530, 
                            445, 395, 440, 370, 405, 315, 185, 170, 140, 145, 75, 40, 40, 35, 45]
                
const AGE_BUCKETS = Dict{Int64, Vector{Int64}}() # key is age group index, value is vector of human indices in that age group
const HOUSEHOLD_MEMBERS = Dict{Int64, Vector{Int64}}() # key is household ID, value is vector of human indices in that household
const FARM_WORKERS = Dict{Int64, Vector{Int64}}() # key is farm ID, value is vector of human indices in that farm

init_state() = init_state(ModelParameters())
function init_state(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    # ip.popsize == 0 && error("no population size given")
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
    # resize the human array to change population size
    resize!(humans, POPSIZE)
    p.flag_init = true # set the flag to true
    #@info "Model initialized with $(POPSIZE) individuals."
    return
end

init_agents() = init_agents(0)
function init_agents(simid) 
    if !p.flag_init
        error("Model is not initialized. Call init_state() first.")
    end
    if p.flag_agents
        @info "Agents already initialized. Resetting households, farming..."
        p.flag_households = false
        p.flag_farming = false
    end
    sum(_POPDISTR) != POPSIZE && error("Population size does not match the sum of the population distribution")
    
    # there is no shuffling going on here, so at index 32250 adults start 
    # this is hardcoded in the model, so don't change it 
    _age_vec = inverse_rle(0:100, _POPDISTR)
    for i in 1:POPSIZE
        humans[i] = Human()
        x = humans[i]
        x.idx = i
        x.age = _age_vec[i]
        x.agegroup = get_age_bucket(_age_vec[i])[1] # get the age group bucket
    end 
   
    # This allocates, but is easier to read and understand
    # the allocations are needed as we need to create the arrays 
    empty!(AGE_BUCKETS) # clear the AGE_BUCKETS dictionary
    for ag in 1:17
        AGE_BUCKETS[ag] = findall(x -> x.agegroup == ag, humans) # filter the humans by age group
    end

    p.flag_agents = true # set the flag to true
    return 
end

function init_households2() 
    if !p.flag_init
        error("Model is not initialized. Call init_state() first.")
    end
    if p.flag_households
        @info "Households already initialized. Run init_agents() first to reinitialize."
        return
    end
    house_ids = create_houseids()

    # we have these IDs separated because the initialized population is not randomized
    # The stateful allows us to iterator over it and take() without resizing arrays
    children_ids = Iterators.Stateful(shuffle!(collect(1:32250))) # children are 0-17 years old, adults are 18+ years old
    adult_ids = Iterators.Stateful(shuffle!(collect(32250:POPSIZE)))
    
    # configuration options for each house of different sizes. Each tuple is (num_adults, num_children)
    comp_options = @SVector [[(1, 0)], 
                    [(1, 1), (2, 0)],
                    [(2, 1), (1, 2), (3, 0)],
                    [(2, 2), (3, 1), (4, 0)],
                    [(2, 3), (3, 2), (4, 1), (5, 0)]]
    ctr = 0
    empty!(HOUSEHOLD_MEMBERS) # clear the HOUSEHOLD_MEMBERS dictionary
    for (i, house_id) in enumerate(house_ids)
        size_code = house_id % 10 # get the size of the house using our encoding scheme
        HOUSEHOLD_MEMBERS[house_id] = Int64[] # initialize the household members dictionary for this house
        all_options = comp_options[size_code] # get the options for this size codecomp_options
        viable = false 
        local o
        while !viable && ctr < 1000000
            ctr += 1
            num_adults, num_child = rand(all_options)
            
            # Check if we have enough adults and children 
            # be careful because e.g., (3, 0) is a valid option, even if there ARE children left or NO children left 
            # (3, 1) is NOT a valid option if there are no children left
            viable = (!isempty(adult_ids) || num_adults == 0) && (!isempty(children_ids) || num_child == 0) 

            if viable
                chunk1 = Iterators.take(children_ids, num_child)
                chunk2 = Iterators.take(adult_ids, num_adults)
                
                # the for loop is really important for materializing the chunks or use collect
                # collect allocates quite a bit of memory, so we use for loop
                for h in chunk1
                    humans[h].hid = house_id
                    push!(HOUSEHOLD_MEMBERS[house_id], h) # add the human index to the household members
                end
                for h in chunk2
                    humans[h].hid = house_id                  
                    push!(HOUSEHOLD_MEMBERS[house_id], h) # add the human index to the household members
                end
            end            
        end
    end

    # Assign adults without a house to a 5+ house randomly. 
    leftover_adults = filter(x -> x.hid == 0, humans)
    large_house_ids = filter(hid -> hid % 10 >= 5, house_ids) # get the house IDs with size 5+
    for h in leftover_adults
        h.hid = rand(large_house_ids) # assign to the last large house
    end

    # error check: no one should be without a house 
    for h in humans
        if h.hid == 0
            @error "Human $(h.idx) has no house assigned!"
        end
    end
    #populate_hh_dict()
    p.flag_households = true # set the flag to true
    return
end 

# for each id, give vector of farm worker IDs
function init_farming() 
    if !p.flag_init
        @error "Model is not initialized. Call init_state() first."
    end
    if p.flag_farming
        @info "Farming already initialized. Run init_agents() first to reinitialize."
        return
    end
    NUM_FARMS = @SVector [131, 85, 28, 12, 3] # number of farms
    EMP_FARMS = @SVector [1:4, 5:9, 10:19, 20:49, 50:99] # number of employees per farm
    total_farms = sum(NUM_FARMS) # total number of farms
    farm_ids = 1:total_farms # farm IDs from 1 to number of farms

    # For each farm, sample the number of employees per farm 
    _employees_per_farm = rand.(EMP_FARMS, NUM_FARMS) # random number of employees per farm
    employees_per_farm = reduce(vcat, _employees_per_farm) # concatenate the arrays into a single array
    total_employees = sum(employees_per_farm) 

    # Now we need to select employees according to an age distribution. 
    # it's possible that doing it this way does not maintain the distribution per farm, but is maintained across all farms. 
    # The relevant age groups are 20-24, 25-34, 35-44, 45-54, 55-64, 65+
    # Sample each employees age group, which gives index from 1 to 6  
    
    PROB_EMP_AG = Categorical(@SVector [0.0375, 0.2375, 0.3250, 0.1875, 0.1250, 0.0875])
    agegroup_of_each_employee = rand(PROB_EMP_AG, total_employees) 
    
    
    #EMP_AG = @SVector [(20, 24), (25, 34), (35, 44), (45, 54), (55, 64), (65, 100)] # age groups
    #indices_per_ag = [shuffle!(collect(get_ag_idx(i[1], i[2]))) for i in EMP_AG] 
    indices_per_ag = [
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 5), # get the indices of humans in each age group
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 6:7),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 8:9),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 10:11),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 12:13),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 14:17),
    ] # get the indices of humans in each age group
    
    # Sample the worker IDs
    sampled_worker_ids = rand.(indices_per_ag[agegroup_of_each_employee])  # ag is array, where each element is an array of indices
    
    # Total farms: 259
    # Total employees needed: 1924
    # Selected unique employees: 1907
    # The reason why we have 1907 employees is because some 
    # employees are selected multiple times due to the random sampling from the age groups.
    # It's too much complexity to ensure that each employee is unique,
    # it doesn't matter because we could've just sampled '1907' total employees to begin with

    # Now we need to assign these sampled workers to farms 
    pos = 1 # position to see who isn't assigned yet from `sampled_worker_ids`
    empty!(FARM_WORKERS) # clear the FARM_WORKERS dictionary
    for (i, fid) in enumerate(farm_ids)
        emp_size =  employees_per_farm[i]
        emp = @view sampled_worker_ids[pos:pos+emp_size-1] # take the first `emp_size` workers from the sampled workers
        for e in emp
            humans[e].fid = fid # assign the farm ID to the human
        end
        pos += emp_size
        FARM_WORKERS[fid] = emp # assign the workers to the farm ID       
    end
    p.flag_farming = true # set the flag to true
    return
end

function get_daily_contacts_farmer(id)
    !(p.flag_households && p.flag_farming) && error("Households and farming must be initialized before getting contacts.")
    farm_id = humans[id].fid 
    farm_id == 0 && error("Human $(id) is not a farmer (no farm ID assigned).")
    house_id = humans[id].hid # get the household ID of the individual
    age_group = humans[id].agegroup # get the age group of the individual
    
    # farm contacts
    all_workers = FARM_WORKERS[farm_id] # get the farm workers for this individual
    num_contact_farm = min(length(all_workers), rand(F1[age_group])) # sample number contacts with farm workers 
    sampled_workers = sample(all_workers, num_contact_farm; replace=false) 
    
    # community contacts
    num_contact_community = rand(F2[age_group]) # sample number of contacts with community individuals
    contact_distr = Categorical(P2[age_group, :]) # distribution of contacts
    sampled_ag_grps = rand(contact_distr, num_contact_community) # sample the age groups for contacts
    sampled_individuals = [rand(AGE_BUCKETS[k]) for k in sampled_ag_grps] # select random individuals from that age buckets
    
    # household members
    sampled_household = HOUSEHOLD_MEMBERS[house_id]
        
    return reduce(vcat, (sampled_workers, sampled_individuals, sampled_household)) # concatenate the arrays of contacts
end

function get_daily_contacts_comm(id)
    house_id = humans[id].hid # get the household ID of the individual
    age_group = humans[id].agegroup # get the age group of the individual
    
    # community contacts
    num_contact_community = rand(F3[age_group]) # sample number of contacts with community individuals
    contact_distr = Categorical(P3[age_group, :]) # distribution of contacts
    sampled_ag_grps = rand(contact_distr, num_contact_community) # sample the age groups for contacts
    sampled_individuals = [rand(AGE_BUCKETS[k]) for k in sampled_ag_grps] # select random individuals from that age buckets
    
    # household members
    sampled_household = HOUSEHOLD_MEMBERS[house_id]
    
    return reduce(vcat, (sampled_individuals, sampled_household)) # concatenate the arrays of contacts
end