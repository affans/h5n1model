function init_model()
    init_state()
    init_agents()
    init_households2()
    init_farming()
    @debug "Model initialized."
end

function create_houseids()
    # creates a unique house ID for each household using a specific encoding scheme
    # where the house id encodes the household size in the last digit
    # formula: houseid = baseid * 10 + size where
    #   base_id is a unique identifier for the house (e.g., sequential number).
    #   size_code is a number from 1 to 5 indicating household size.
    # to find household size from a house_id: size = house_id % 10
    # to get the unique baseid: base_id = house_id // 10

    # the houses here can accomodate 142285 individuals, we have a fixed size of POPSIZE = 153525,
    # so we have less houses than popsize - but the last house size is 5+ so throw the rest in there
    household_sizes = @SVector [13180, 16335, 7265, 7660, 8800] #sum 53240
    sizes = inverse_rle(1:5, household_sizes) # size of each house
    bases = 1:sum(household_sizes) # unique ID of each house
    houseids = bases .* 10 .+ sizes
    return houseids
end

function poplast!(a, n=1)
    # pops the last n elements from an array a,
    # resizing the array in place and returning the popped elements
    # initially used but stopped since it's allocating -- Use Iterators.take instead
    n = min(n, length(a))
    popped = a[end-n+1:end]
    resize!(a, length(a) - n)
    return popped
end

function get_age_bucket(x)
    # returns the age bucket index (from 1 to 17) and the range for a given age x
    # the age buckets are defined based on google sheets for contacts
    buckets = @SVector([(0, 4), (5, 9), (10, 14), (15, 19), (20, 24), (25, 29),
                        (30, 34), (35, 39), (40, 44), (45, 49), (50, 54), (55, 59),
                        (60, 64), (65, 69), (70, 74), (75, 79), (80, 100)])
    if x >= 80
        bucket_idx = 17  # last bucket index
    else
        bucket_idx = div(x, 5) + 1  # +1 for 1-based index
        range_start = (bucket_idx - 1) * 5
        range_end = range_start + 4
    end
    return bucket_idx, buckets[bucket_idx]
end

# Contacts Distributions
####################

# Distribution for number of contacts for FARMER at WORK
# The index of each element represents the age group (from 1 to 17 - see `get_age_bucket`)
const F1 = @SVector [
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(1.187694257, 0.1633981111),
    NegativeBinomial(1.394465841, 0.1586527642),
    NegativeBinomial(1.558780229, 0.1748411908),
    NegativeBinomial(1.430978535, 0.1742133948),
    NegativeBinomial(1.170066711, 0.1389676949),
    NegativeBinomial(1.165204064, 0.1389676949),
    NegativeBinomial(1.054874243, 0.1439876496),
    NegativeBinomial(0.9672436716, 0.1718948873),
    NegativeBinomial(0.6923976042, 0.232822678),
    NegativeBinomial(0.3955447194, 0.2554829004),
    NegativeBinomial(0.1134947043, 0.3225621743),
    NegativeBinomial(0.2443440758, 0.3365528229),
    NegativeBinomial(0.02440799928, 0.3816104412)
]


# Distribution for number of contacts for FARMER with COMMUNITY
# The index of each element represents the age group (from 1 to 17 - see `get_age_bucket`)
const F2 = @SVector [
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.0100000000, 0.06317366633), ## will not be used
    NegativeBinomial(0.3359411191, 0.1633981111),
    NegativeBinomial(0.2244125912, 0.1586527642),
    NegativeBinomial(0.2887244332, 0.1748411908),
    NegativeBinomial(0.3248316133, 0.1742133948),
    NegativeBinomial(0.2689281652, 0.1389676949),
    NegativeBinomial(0.2314521561, 0.1389676949),
    NegativeBinomial(0.2194465676, 0.1439876496),
    NegativeBinomial(0.4066662282, 0.1718948873),
    NegativeBinomial(0.5251787696, 0.232822678),
    NegativeBinomial(0.8122052032, 0.2554829004),
    NegativeBinomial(1.101686717, 0.3225621743),
    NegativeBinomial(1.046743663, 0.3365528229),
    NegativeBinomial(1.180641593, 0.3816104412)
]




# Distribution for number of contacts for General
# The index of each element represents the age group (from 1 to 17 - see `get_age_bucket`)
const F3 = @SVector [
    NegativeBinomial(0.99406, 0.25556), 
    NegativeBinomial(1.50982, 0.23675),
    NegativeBinomial(1.75702, 0.21788),
    NegativeBinomial(1.84065, 0.22767),
    NegativeBinomial(1.50252, 0.16340),
    NegativeBinomial(1.53524, 0.15865),
    NegativeBinomial(1.87315, 0.17484),
    NegativeBinomial(1.80414, 0.17421),
    NegativeBinomial(1.50041, 0.13897),
    NegativeBinomial(1.53205, 0.13897),
    NegativeBinomial(1.37195, 0.14399),
    NegativeBinomial(1.34388, 0.17189),
    NegativeBinomial(1.18221, 0.23282),
    NegativeBinomial(1.14081, 0.25548),
    NegativeBinomial(1.19670, 0.32256),
    NegativeBinomial(1.29392, 0.33655),
    NegativeBinomial(1.20159, 0.38161)
]

# Distribution of Who meets Who for FARMER at COMMUNITY
# The index of each element represents the age group (from 1 to 17 - see `get_age_bucket`)
const P2 = @SVector [
    Categorical([1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), ## will not be used
    Categorical([1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), ## will not be used
    Categorical([1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), ## will not be used
    Categorical([1.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), ## will not be used
    Categorical([0.0006, 0.0028, 0.0327, 0.0994, 0.2464, 0.1490, 0.0970, 0.0540, 0.0788, 0.0781, 0.0403, 0.0436, 0.0298, 0.0231, 0.0079, 0.0131, 0.0034]),
    Categorical([0.0073, 0.0013, 0.0098, 0.0296, 0.1622, 0.2926, 0.1495, 0.0931, 0.0524, 0.0318, 0.0529, 0.0477, 0.0334, 0.0170, 0.0059, 0.0127, 0.0008]),
    Categorical([0.0294, 0.0066, 0.0207, 0.0380, 0.0591, 0.1763, 0.1785, 0.0988, 0.0718, 0.0893, 0.0709, 0.0798, 0.0346, 0.0270, 0.0136, 0.0025, 0.0031]),
    Categorical([0.0275, 0.0202, 0.0269, 0.0433, 0.0720, 0.1162, 0.1606, 0.1593, 0.1033, 0.0758, 0.0562, 0.0446, 0.0296, 0.0254, 0.0246, 0.0067, 0.0078]),
    Categorical([0.0095, 0.0284, 0.0256, 0.0424, 0.0855, 0.1215, 0.1275, 0.1048, 0.1189, 0.0985, 0.0766, 0.0702, 0.0297, 0.0242, 0.0184, 0.0108, 0.0075]),
    Categorical([0.0098, 0.0043, 0.0918, 0.0282, 0.0839, 0.0657, 0.0718, 0.0779, 0.0934, 0.0986, 0.1104, 0.1009, 0.0510, 0.0442, 0.0385, 0.0226, 0.0070]),
    Categorical([0.0069, 0.0287, 0.0144, 0.0508, 0.0826, 0.0909, 0.0736, 0.0850, 0.0641, 0.0879, 0.0926, 0.1178, 0.0734, 0.0487, 0.0301, 0.0268, 0.0257]),
    Categorical([0.0166, 0.0108, 0.0105, 0.0358, 0.0656, 0.0985, 0.0891, 0.0659, 0.0761, 0.0675, 0.1130, 0.1197, 0.1001, 0.0599, 0.0309, 0.0173, 0.0227]),
    Categorical([0.0090, 0.0168, 0.0183, 0.0259, 0.0658, 0.0785, 0.0883, 0.0630, 0.0577, 0.0526, 0.0895, 0.1105, 0.1204, 0.0894, 0.0629, 0.0217, 0.0297]),
    Categorical([0.0076, 0.0047, 0.0054, 0.0145, 0.0398, 0.0554, 0.0771, 0.0522, 0.0507, 0.0619, 0.0936, 0.1137, 0.0950, 0.1347, 0.1065, 0.0430, 0.0442]),
    Categorical([0.0032, 0.0047, 0.0049, 0.0163, 0.0430, 0.0549, 0.0788, 0.0682, 0.0419, 0.0700, 0.0917, 0.1068, 0.0928, 0.1100, 0.1011, 0.0521, 0.0596]),
    Categorical([0.0022, 0.0064, 0.0039, 0.0134, 0.0350, 0.0615, 0.0974, 0.0470, 0.0454, 0.0773, 0.0788, 0.0797, 0.0649, 0.0940, 0.0885, 0.1003, 0.1043]),
    Categorical([0.0032, 0.0022, 0.0029, 0.0210, 0.0083, 0.0530, 0.0525, 0.0723, 0.0405, 0.0479, 0.1196, 0.0804, 0.0828, 0.1048, 0.0815, 0.0932, 0.1339])
]


# Distribution of Who meets Who for NON-FARMER at COMMUNITY
# The index of each element represents the age group (from 1 to 17 - see `get_age_bucket`)
const P3 = @SVector [
    Categorical([0.443, 0.097, 0.015, 0.011, 0.029, 0.044, 0.086, 0.060, 0.060, 0.041, 0.042, 0.033, 0.017, 0.011, 0.006, 0.002, 0.003]),
    Categorical([0.080, 0.468, 0.070, 0.011, 0.012, 0.063, 0.043, 0.061, 0.054, 0.038, 0.066, 0.010, 0.010, 0.006, 0.005, 0.002, 0.001]),
    Categorical([0.005, 0.082, 0.504, 0.047, 0.017, 0.047, 0.114, 0.048, 0.043, 0.034, 0.022, 0.012, 0.014, 0.005, 0.004, 0.001, 0.001]),
    Categorical([0.006, 0.005, 0.133, 0.215, 0.082, 0.072, 0.137, 0.101, 0.067, 0.054, 0.040, 0.041, 0.013, 0.009, 0.007, 0.012, 0.006]),
    Categorical([0.006, 0.010, 0.022, 0.053, 0.177, 0.128, 0.139, 0.099, 0.097, 0.105, 0.052, 0.044, 0.026, 0.019, 0.011, 0.007, 0.005]),
    Categorical([0.007, 0.028, 0.028, 0.028, 0.101, 0.147, 0.143, 0.126, 0.098, 0.073, 0.073, 0.066, 0.031, 0.023, 0.011, 0.007, 0.010]),
    Categorical([0.010, 0.029, 0.065, 0.049, 0.086, 0.109, 0.121, 0.112, 0.094, 0.077, 0.074, 0.075, 0.035, 0.029, 0.016, 0.011, 0.008]),
    Categorical([0.015, 0.028, 0.049, 0.034, 0.064, 0.105, 0.129, 0.142, 0.118, 0.087, 0.077, 0.077, 0.027, 0.020, 0.014, 0.005, 0.009]),
    Categorical([0.027, 0.065, 0.047, 0.029, 0.076, 0.093, 0.116, 0.124, 0.099, 0.093, 0.081, 0.078, 0.029, 0.018, 0.010, 0.006, 0.009]),
    Categorical([0.045, 0.030, 0.056, 0.039, 0.090, 0.078, 0.095, 0.106, 0.103, 0.088, 0.089, 0.085, 0.035, 0.027, 0.017, 0.010, 0.007]),
    Categorical([0.027, 0.031, 0.036, 0.034, 0.070, 0.074, 0.097, 0.094, 0.094, 0.091, 0.091, 0.096, 0.057, 0.043, 0.027, 0.017, 0.021]),
    Categorical([0.019, 0.022, 0.019, 0.027, 0.054, 0.086, 0.100, 0.098, 0.086, 0.092, 0.105, 0.107, 0.068, 0.048, 0.030, 0.018, 0.021]),
    Categorical([0.019, 0.015, 0.020, 0.022, 0.052, 0.081, 0.106, 0.081, 0.064, 0.069, 0.105, 0.115, 0.095, 0.064, 0.046, 0.020, 0.026]),
    Categorical([0.005, 0.003, 0.004, 0.013, 0.046, 0.070, 0.099, 0.073, 0.065, 0.065, 0.092, 0.109, 0.082, 0.110, 0.085, 0.037, 0.042]),
    Categorical([0.003, 0.006, 0.005, 0.016, 0.043, 0.055, 0.078, 0.069, 0.045, 0.070, 0.090, 0.106, 0.094, 0.108, 0.102, 0.051, 0.059]),
    Categorical([0.002, 0.005, 0.031, 0.025, 0.033, 0.055, 0.085, 0.043, 0.050, 0.076, 0.079, 0.081, 0.073, 0.095, 0.083, 0.091, 0.093]),
    Categorical([0.003, 0.002, 0.003, 0.021, 0.008, 0.052, 0.052, 0.071, 0.042, 0.047, 0.117, 0.085, 0.084, 0.103, 0.083, 0.092, 0.135])
]