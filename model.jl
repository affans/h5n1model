using StatsBase
using Random
using Distributions
using DelimitedFiles
using StaticArrays

# make sure the order doesn't change
@enum INFTYPE UNDEF = 0 SUS = 1 EXP = 2 ASYMP = 3 SYMP = 4 REC = 5 # infection status
@enum VACSCN1 NONE FARMONLY FARM_AND_HH 
@enum VACSCN2 A0 A1 A2 A3 A4 # vaccination scenarios, # A0 = No scenario 

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
    iso::Bool = false # isolation flag
    vaceff1::Float32 = 0.0 
    vaceff2::Float32 = 0.0
    totalinfect::Int64 = 0 # total number of infections this individual has had
    timeofinfect::Int64 = 0 # time of infection
    # eff::Float64 = 0.0 # vaccine efficacy
    # exp::Int64 = 0 # waning time of vaccine
end
Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## system parameters
# Base.@kwdef mutable struct ModelParameters
#     flag_init = false # flag to indicate if the model is initialized
# end

include("helpers.jl") # include the helper functions
include("debug.jl")

# constant variables
const POPSIZE = 153565
const humans = Array{Human}(undef, POPSIZE)
#const p = ModelParameters()  ## setup default parameters
const AGE_BUCKETS = Dict{Int64,Vector{Int64}}() # key is age group index, value is vector of human indices in that age group
const HOUSEHOLD_MEMBERS = Dict{Int64,Vector{Int64}}() # key is household ID, value is vector of human indices in that household
const FARM_WORKERS = Dict{Int64,Vector{Int64}}() # key is farm ID, value is vector of human indices in that farm
const INF_PERIOD = Gamma(3.2949, 2.7314)
const EXP_PERIOD = Weibull(2.3015, 3.7242) # incubation period
const will_meet = falses(POPSIZE) # array to mark who will meet whom

function init_state()
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    # ip.popsize == 0 && error("no population size given")
    # for x in propertynames(p)
    #     setfield!(p, x, getfield(ip, x))
    # end
    # resize the human array to change population size
    #resize!(humans, POPSIZE)
    return
end

init_agents() = init_agents(0)
function init_agents(simid)
    _POPDISTR = @SVector [1530, 1545, 1565, 1670, 1700, 1820, 1755, 1790, 1820, 1935, 1875, 1925,
        1925, 1910, 1890, 1800, 1895, 1900, 1820, 1855, 1995, 2135, 2040, 2060,
        2030, 2120, 1985, 2020, 2125, 2140, 2120, 2230, 2240, 2065, 2105, 2190,
        2160, 2195, 2225, 2210, 2190, 2025, 2005, 1965, 1845, 1870, 1920, 1840,
        1835, 1855, 1880, 1895, 1870, 1915, 1645, 1820, 1885, 1970, 1875, 1755,
        1850, 1945, 1845, 1850, 1790, 1720, 1810, 1735, 1595, 1520, 1460, 1515,
        1430, 1500, 1410, 1135, 1060, 1050, 925, 860, 820, 735, 730, 755, 570, 530,
        445, 395, 440, 370, 405, 315, 185, 170, 140, 145, 75, 40, 40, 35, 45]

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

    # empty the dicts 
    empty!(AGE_BUCKETS) # clear the AGE_BUCKETS dictionary
    empty!(HOUSEHOLD_MEMBERS) # clear the HOUSEHOLD_MEMBERS dictionary
    empty!(FARM_WORKERS) # clear the FARM_WORKERS dictionary

    populate_age_dict() 
    return
end

function populate_age_dict() 
    # This allocates, but is easier to read and understand
    # the allocations are needed as we need to create the arrays 
    empty!(AGE_BUCKETS) # clear the AGE_BUCKETS dictionary
    for ag in 1:17
        AGE_BUCKETS[ag] = findall(x -> x.agegroup == ag, humans) # filter the humans by age group
    end
end

function init_households2()
    # if !p.flag_init
    #     error("Model is not initialized. Call init_state() first.")
    # end

    house_ids = create_houseids()
    
    # we have these IDs separated because the initialized population is not randomized
    # The stateful allows us to iterator over it and take() without resizing arrays
    children_ids = Iterators.Stateful(shuffle!(collect(1:32250))) # children are 0-17 years old, adults are 18+ years old
    adult_ids = Iterators.Stateful(shuffle!(collect(32250:POPSIZE)))

    #children_ids = shuffle!(collect(1:32250)) # children are 0-17 years old, adults are 18+ years old
    #adult_ids = shuffle!(collect(32250:POPSIZE))
    # configuration options for each house of different sizes. Each tuple is (num_adults, num_children)
    comp_options = [
        [(1, 0)],
        [(1, 1), (2, 0)],
        [(2, 1), (1, 2), (3, 0)],
        [(2, 2), (3, 1), (4, 0)],
        [(2, 3), (3, 2), (4, 1), (5, 0)]
    ]
    ctr = 0
    
    for (i, house_id) in Iterators.enumerate(house_ids)
        #HOUSEHOLD_MEMBERS[house_id] .= 0 # initialize the household members dictionary for this house
        size_code = house_id % 10 # get the size of the house using our encoding scheme  
        all_options = comp_options[size_code] # get the options for this size codecomp_options
        viable = false
        while !viable && ctr < 1000000
            ctr += 1
            num_adults, num_child = rand(all_options)
            # Check if we have enough adults and children 
            # be careful because e.g., (3, 0) is a valid option, even if there ARE children left or NO children left 
            # (3, 1) is NOT a valid option if there are no children left
            # NOTE: viable does not mean there are exactly that many adults and children,
            # it means that there are enough adults and children to partially fill the house
            # so if viable (3, 2) and there are 3 adults and 1 child left, it is still viable 
            # it's rare, but often the case we run out of children
            viable = (!isempty(adult_ids) || num_adults == 0) && (!isempty(children_ids) || num_child == 0)
           
            if viable
                chunk1 = Iterators.take(children_ids, num_child)
                chunk2 = Iterators.take(adult_ids, num_adults)
                
                # the for loop is really important for materializing the chunks or use collect
                # collect allocates quite a bit of memory, so we use for loop
                for h in chunk1
                    humans[h].hid = house_id
                end
                for h in chunk2
                    humans[h].hid = house_id
                end
            end
        end
    end

    # Assign adults without a house to a 5+ house randomly. 
    leftover_adults = filter(x -> x.hid == 0, humans)
    large_house_ids = filter(hid -> hid % 10 >= 5, house_ids) # get the house IDs with size 5+
    for h in leftover_adults
        h.hid = rand(large_house_ids) # assign to the last large house
        #push!(HOUSEHOLD_MEMBERS[house_id], h)
    end

    #error check: no one should be without a house 
    for h in humans
        if h.hid == 0
            # @error "Human $(h.idx) has no house assigned!"
            #return
        end
    end

    # populate the HOUSEHOLD_MEMBERS dictionary
    populate_household_dict()

    return
end

function populate_household_dict() 
    houseidx = [h.hid for h in humans] # get the household IDs of all humans
    empty!(HOUSEHOLD_MEMBERS) # clear the HOUSEHOLD_MEMBERS dictionary
    for k in unique(houseidx) # for each unique household ID
        HOUSEHOLD_MEMBERS[k] = Int64[]
    end
    for (i, v) in enumerate(houseidx)
        push!(HOUSEHOLD_MEMBERS[v], i) # add the human index to the household members
    end
    return
end

function init_farming()
    NUM_FARMS = @SVector [85, 28, 12, 3, 3] # number of farms
    #NUM_FARMS = @SVector [131, 28, 12, 3, 3] # number of farms
    #@info "Total number of farms: $(sum(NUM_FARMS))"
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

    indices_per_ag = [
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 5), # get the indices of humans in each age group
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 6:7),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 8:9),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 10:11),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 12:13),
        mapfoldl(i -> AGE_BUCKETS[i], vcat, 14:17),
    ] # get the indices of humans in each age group

    how_many_each_ag = countmap(agegroup_of_each_employee)
    sampled_worker_ids = zeros(Int64, total_employees) # array to hold the sampled worker IDs

    a1 = sample(indices_per_ag[1], how_many_each_ag[1], replace=false) 
    a2 = sample(indices_per_ag[2], how_many_each_ag[2], replace=false)
    a3 = sample(indices_per_ag[3], how_many_each_ag[3], replace=false)
    a4 = sample(indices_per_ag[4], how_many_each_ag[4], replace=false)
    a5 = sample(indices_per_ag[5], how_many_each_ag[5], replace=false)
    a6 = sample(indices_per_ag[6], how_many_each_ag[6], replace=false)
    sampled_worker_ids = reduce(vcat, [a1, a2, a3, a4, a5, a6]) # concatenate the arrays into a single array
    #sampled_worker_ids = rand.(indices_per_ag[agegroup_of_each_employee])  # ag is array, where each element is an array of indices

    # Now we need to assign these sampled workers to farms 
    pos = 1 # position to see who isn't assigned yet from `sampled_worker_ids`
    empty!(FARM_WORKERS) # clear the FARM_WORKERS dictionary
    for (i, fid) in enumerate(farm_ids)
        emp_size = employees_per_farm[i]
        emp = @view sampled_worker_ids[pos:pos+emp_size-1] # take the first `emp_size` workers from the sampled workers
        for e in emp
            humans[e].fid = fid # assign the farm ID to the human
        end
        pos += emp_size
        FARM_WORKERS[fid] = emp # assign the workers to the farm ID       
    end
    return
end

# function insert_infection()
#     farm = rand(keys(FARM_WORKERS))
#     idx = rand(FARM_WORKERS[farm]) # sample a random farm worker
#     x = humans[idx] # get the human object
#     x.timeofinfect = 1
#     x.swap = SYMP
#     activate_swaps(x) # activate the swap to symptomatic
#     return idx # return the indices of infected individuals
# end
insert_infection() = insert_infection(1)
function insert_infection(num_of_infections)
    farm = rand(keys(FARM_WORKERS), num_of_infections)
    ids = zeros(Int64, num_of_infections)
    for (i, fid) in enumerate(farm)
        idx = rand(FARM_WORKERS[fid]) # sample a random farm worker
        x = humans[idx] # get the human object
        x.timeofinfect = 1
        x.swap = SYMP
        activate_swaps(x) # activate the swap to symptomatic
        ids[i] = idx
    end
    return ids
end

### TIME STEP FUNCTIONS
function time_loop(;simid=1, 
                    beta=0.0,
                    init_inf=0,
                    iso_day=-1, 
                    iso_prop=0.0, 
                    vac_scn1=NONE, 
                    vac_scn2=A0,
                    vac_cov=0.0
                    )
    Random.seed!(simid*293)
    init_model()
    insert_infection(init_inf) # insert an infection in the farm
    
    # data collection variables
    max_time = 365
    incidence_hh = zeros(Int64, max_time) # incidence per day in household
    incidence_fm = zeros(Int64, max_time) # incidence per day in farm
    incidence_cm = zeros(Int64, max_time) # incidence per day in comm

    for t in Base.OneTo(max_time)
        for (i, x) in enumerate(humans)
            iso_dynamics(x, iso_day, iso_prop) # check if the individual needs to be isolated or not (before natural history) 
            natural_history(x) # move through the natural history of the disease first         
            _, tih, tif, tic = transmission_with_contacts!(simid, t, x, beta)
            activate_swaps(x) # essentially "midnight"
            
            # data collection 
            incidence_cm[t] += tic # add the community incidence
            incidence_fm[t] += tif # add the farm incidence
            incidence_hh[t] += tih # add the household incidence
        end
        
        # vaccine init at 6 weeks later
        if t == 42 && vac_scn1 != NONE # if it's the first month and we have a vaccination scenario 
            @debug "Day $(t) of simulation $(simid)."
            @debug "Vaccination scenario: $(vac_scn1), $(vac_scn2) with coverage: $(vac_cov)."
            init_vac(simid, vac_scn1, vac_scn2, vac_cov)
        end
    end

    # who infect who vector for  R-effective calculation.
    # who_infect_who = [(x.totalinfect, x.timeofinfect) for x in humans if x.inf == REC]
    return incidence_hh, incidence_fm, incidence_cm #, who_infect_who
end

function calibrate(beta; iso_day=-1, iso_prop=0.0) 
    # initialize the model
    init_model()
    init_infect_id = insert_infection(1)[1] # get the first infected individual
    x = humans[init_infect_id] # get the infected individual
    max_inf_time = x.st + 1
    total_infected = 0 
    for t in 1:max_inf_time
        iso_dynamics(x, iso_day, iso_prop) # check if the individual needs to be isolated or not (before natural history) 
        natural_history(x) # move through the natural history of the disease first
        _, tih, tif, tic = transmission_with_contacts!(1, t, x, beta)
        total_infected += (tih + tif + tic) # count the total infected individuals
        @debug "End Day $(t): tis: $(x.tis), st: $(x.st), swap: $(x.swap)  met $(tm), inf total: $(total_infected), inf hh $(tih), inf farm $(tif), inf comm $(tic)."
        activate_swaps(x)
    end
    @debug "Time end: $(x.idx), tis: $(x.tis), st: $(x.st), swap: $(x.swap), st: $(x.st)"
    @debug "Total infected: $(total_infected) individuals."
    return total_infected
end

function natural_history(x::Human)
    x.tis += 1
    will_swap = false
    if x.tis == x.st # move through nat history of disease
        will_swap = true
        x.inf ∈ (SUS, REC) && error("SUS/REC state can not expire")
        if x.inf == EXP 
            prob_of_symp = 0.97 * (1.0 - x.vaceff2)
            if rand() < prob_of_symp
                x.swap = SYMP # swap to symptomatic
            else
                x.swap = ASYMP # 3% of incubating individuals become asymptomatic
            end
        end
        if x.inf ∈ (SYMP, ASYMP) # symptomatic or asymptomatic
            x.swap = REC
        end
    end
    return will_swap
end

function is_infectious(x::Human)
    return x.inf ∈ (SYMP, ASYMP)
end

@inline function iso_dynamics(x::Human, iso_day, iso_prop)
    (x.inf == SYMP && x.tis == iso_day && rand() < iso_prop) && (x.iso = true)
end

function activate_swaps(x::Human)
    if x.swap ≠ UNDEF
        x.tis = 0  # reset time in state
        if x.swap == EXP
            x.inf = EXP # swap to incubating
            x.st = round(Int64, rand(EXP_PERIOD)) # set the incubation period
        elseif x.swap == ASYMP
            x.inf = ASYMP # swap to asymptomatic
            inf_per = max(round(Int64, rand(INF_PERIOD)), 0)
            x.st = round(Int64, inf_per) # set the infection period
        elseif x.swap == SYMP
            x.inf = SYMP # swap to symptomatic
            inf_per = max(round(Int64, rand(INF_PERIOD)), 0)
            x.st = round(Int64, inf_per) # set the infection period
        elseif x.swap == REC
            x.inf = REC # swap to recovered
            x.st = typemax(Int16) # set the recovery period to max
        end
        x.swap = UNDEF # reset the swap state
    end
    return
end

function init_vac(simid, scn1::VACSCN1, scn2::VACSCN2, coverage) 
    # at 4 weeks from initial infection, 
    # vaccinate (i) farm workers with some coverage, including their family members 
    rng2 = MersenneTwister(simid + 2638)
    if scn2 == A0 
        eff1 = 0.0 # no protection against infection
        eff2 = 0.0 # no additional protection for symptomatic disease
    elseif scn2 == A1 
        eff1 = 0.50 # 50% protection against infection
        eff2 = 0.0 # no additional protection for symptomatic disease
    elseif scn2 == A2
        eff1 = 0.0 # 50% protection against infection
        eff2 = 0.85 # 50% additional protection for symptomatic disease
    elseif scn2 == A3
        eff1 = 0.50 # no protection against infection
        eff2 = 0.85 # 50% additional protection for symptomatic disease
    end
    
    for (k, v) in FARM_WORKERS
        if rand(rng2) < coverage
            for idx in v
                x = humans[idx]
                if scn1 == FARMONLY 
                    x.vaceff1 = eff1 # mark the individual as vaccinated
                    x.vaceff2 = eff2
                elseif scn1 == FARM_AND_HH
                    hh = HOUSEHOLD_MEMBERS[x.hid] # get the household members of the farm worker 
                    for hidx in hh
                        y = humans[hidx] # get the household member
                        if y.age >= 1
                            y.vaceff1 = eff1 # mark the individual as vaccinated
                            y.vaceff2 = eff2
                        end
                    end
                end
            end 
        end
    end
    # since this uses random numbers in the middle of the simulation, we need to reset the random number generator
    return
end

@enum BETA_REDUC_TYPE BR_HH BR_FARM BR_COMM # beta reduction types
function transmission_with_contacts!(simid, time, x::Human, beta)
    # This is an allocation free method to get daily contacts for an individual
    # for each contact, we check for transmission right away to avoid allocating 
    # however, the flag on the global will_meet will set to true for each contact
    # (which would reset for the next infected individual)
    !(x.inf ∈ (ASYMP, SYMP)) && return (0, 0, 0, 0) # if the individual is not infected, return 0 contacts and 0 infections
    
    farm_id = x.fid
    house_id = x.hid # get the household ID of the individual
    age_group = x.agegroup # get the age group of the individual
    
    # counting variables
    total_inf_household = 0
    total_inf_farms = 0
    total_inf_community = 0 
    total_meet = 0 

    # not sure if this is used anymore -- defined at the global level
    # will_meet .= false # reset the will meet array for this individual

    # household contacts
    for idx in HOUSEHOLD_MEMBERS[house_id]
        # will_meet[idx] = true # mark the household member as will meet (in case needed later)
        total_inf_household += check_for_transmission(time, x, humans[idx], beta, BR_HH) # check for transmission
        total_meet += 1
    end
    
    # if the individual is isolated, return
    if x.iso 
        return total_meet, total_inf_household, total_inf_farms, total_inf_community
    end

    # farm contacts
    if farm_id > 0 # if the individual is a farm worker and not isolated
        all_workers = shuffle!(FARM_WORKERS[farm_id]) # get the farm workers for this individual
        num_contact_farm = min(length(all_workers), rand(F1[age_group])) # sample number contacts with farm workers 
        sampled_workers = @view all_workers[1:num_contact_farm] # take the first num_contact_farm workers since they are shuffled
        #@debug "Farmer $(x.idx) will meet $(num_contact_farm) farm workers."
        for idx in sampled_workers
            #@debug "Farmer $(x.idx) will meet farm worker $(idx)." 
            # will_meet[idx] = true # mark the worker as will meet
            total_inf_farms += check_for_transmission(time, x, humans[idx], beta, BR_FARM) # check for transmission
            total_meet += 1
        end
    end

    # community contacts
    if farm_id > 0
        num_contact_community = rand(F2[age_group]) # sample number of contacts with community individuals
        contact_distr = P2[age_group] # distribution of contacts, allocates
    else
        num_contact_community = rand(F3[age_group]) # sample number of contacts with community individuals
        contact_distr = P3[age_group] # distribution of contacts
    end
    for i in Base.OneTo(num_contact_community) # go through each contact
        k = rand(contact_distr) # sample the age group for contacts
        idx = rand(AGE_BUCKETS[k]) # sample a random person from the age bucket for that age group
        # will_meet[idx] = true # mark the individual in the age group as will meet
        total_inf_community += check_for_transmission(time, x, humans[idx], beta, BR_COMM) # check for transmission
        total_meet += 1
        #@debug "Farmer $(x.idx) will meet community individual $(idx) from age group $(k)."
    end
  
    return total_meet, total_inf_household, total_inf_farms, total_inf_community # return the total number of contacts and infections
end

function check_for_transmission(time, x::Human, y::Human, beta, br::BETA_REDUC_TYPE)
    infect = 0
    effbeta = beta
    if br == BR_HH 
        effbeta *= 0.43 * beta
    elseif br == BR_FARM 
        effbeta *= 0.23 * beta
    elseif br == BR_COMM # community transmission
        effbeta *= 0.20 * beta
    end

    if x.inf == ASYMP # if the individual is symptomatic 
        effbeta = effbeta * 0.20
    end

    effbeta = effbeta * (1.0 - y.vaceff1) # vaccine efficacy against infection

    if y.inf == SUS && y.swap == UNDEF # if the individual is susceptible
        if rand() < effbeta # check if transmission occurs
            x.totalinfect += 1 
            y.timeofinfect = time
            y.swap = EXP # swap to incubating
            infect = 1
        end
    end
    return infect
end
