function poplast!(a, n=1)
    n = min(n, length(a))
    popped = a[end-n+1:end]
    resize!(a, length(a) - n)
    return popped
end

# get's the array indices for a given age or age range
# this is faster/nonallocating than using `findall` or `filter`
@inline function get_ag_idx(n1::Int64, n2=n1) 
    # returns the age group index for a given age
    if n1 < 0 || n1 > 100 || n2 < 0 || n2 > 100
        error("Age $n1, $n2 is out of bounds (0-100)")
    end
    pop = cumsum(_POPDISTR) # cumulative sum of the population distribution
    start_idx = n1 == 0 ? 1 : (pop[n1] + 1)
    end_idx = pop[n2+1]
    return start_idx:end_idx 
end