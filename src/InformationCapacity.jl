module InformationCapacity
export xcov, autocov, chunkwise, legendre_prods, enumerate_polys_bydeg, enumerate_polys,
    capacity, x2X, sample, find_threshold, eval_capacities, analyze_capacities,
    optimize_offset, bisect_int, getYs,
    DOPO
    
using Jacobi


include("capacity.jl")
include("dopo.jl")

end # module
