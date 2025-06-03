function poplast!(a, n=1)
    n = min(n, length(a))
    popped = a[end-n+1:end]
    resize!(a, length(a) - n)
    return popped
end