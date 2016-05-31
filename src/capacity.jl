
using PyPlot
using LaTeXStrings 
plt[:style][:use]("ggplot")


function xcov(m1, m2)
    n1, p1 = size(m1)
    n2, p2 = size(m2)

    p = min(p1,p2)

    YXT = zeros(n1,n2)

    for kk=1:p
        for jj=1:n2
            @simd for ll=1:n1
                @inbounds YXT[ll,jj] += m1[ll,kk] * m2[jj,kk]
            end
        end
    end
    scale!(YXT, 1./p)

    YXT
end

function autocov(m)
    xcov(m,m)
end

immutable Chunked{T}
    seq::Vector{T}
    clen::Int
    cinc::Int
end

import Base:getindex, size, convert

@inline getindex(c::Chunked, k,l) = (k==1 ? 1 : c.seq[k-1+(l-1)*c.cinc])
size(c::Chunked) = (1+c.clen, ceil(Int, (length(c.seq)-c.clen + 1)/c.cinc))
size(c::Chunked, jj::Int) = size(c)[jj]
function convert{T}(::Type{Matrix}, c::Chunked{T})
    ret = zeros(T, size(c)...)
    for kk=1:size(c,2)
        for ll=1:size(c,1)
            ret[ll,kk] = c[ll,kk]
        end
    end
    ret
end

function x2X(x, chunklen, chunkinc=nothing)
    if chunkinc == nothing
        chunkinc=chunklen
    end
    N = length(x) - chunklen + 1
    if chunkinc <= chunklen
        ns = ceil(Int, N/chunkinc)
    else
        ns = floor(Int, N/chunkinc)
    end
    X = ones(chunklen+1, ns)

    for kk=0:ns-1
        X[2:end,kk+1] = x[kk*chunkinc+1:kk*chunkinc+chunklen]
    end
    X
end



function chunkwise(seq, chunklen, chunkinc=nothing)
    if chunkinc == nothing
        chunkinc=chunklen
    end
    n = length(seq)
    Chunked(seq, chunklen, chunkinc)
end

function legendre_prods(u, djks)
    T = length(u)
    depth, nks = size(djks)
    maxdeg = maximum(djks)

    polys = zeros(eltype(u), T, maxdeg)
    for kk=1:maxdeg
        legendre!(u, kk, sub(polys, 1:T, kk))
    end
    res = ones(eltype(u), T-depth+1, nks)
    for kk=1:nks
        for tau=1:depth
            if djks[tau, kk] > 0
                res[:,kk] .*= polys[tau:end-(depth-tau), djks[tau, kk]]
            end
        end
    end
    res'
end

_bydeg_cache = Dict{Tuple{Int,Int},Matrix{Int}}()
function enumerate_polys_bydeg(depth, deg)
    if (depth, deg) in keys(_bydeg_cache)
        return _bydeg_cache[(depth, deg)]
    end
    if depth < 1
        error("Need at least one time slice")
    end
    if depth == 1
        return ones(Int,1,1)*deg
    end
    if deg == 0
        return zeros(Int, depth, 1)
    end
    res_by_l = []
    for l=0:deg
        rec = enumerate_polys_bydeg(depth-1, deg-l)
        rec = [rec; l*ones(Int, 1, size(rec, 2))]
        push!(res_by_l, rec)
    end
    res = hcat(res_by_l...)
    _bydeg_cache[(depth, deg)] = res
    return res

end


function enumerate_polys(depth, maxdeg)
    ret = zeros(Int, depth, 0)
    for d=1:maxdeg
        ret = hcat(ret, enumerate_polys_bydeg(depth, d))
    end
    ret
end

function estimate_kernel(Y, Xs, XcorrFact=nothing)
    if XcorrFact === nothing
        XcorrFact = cholfact(autocov(Xs))
    end
    XYt = xcov(Xs,Y)
    kernel = (XcorrFact \ XYt)'
end


function capacity(Y::Matrix, Xs, XcorrFact=nothing)
    if XcorrFact === nothing
        Xac = autocov(Xs)
        XcorrFact = lufact(Xac)
    end
#     println(size(Xs))
#     println(size(Y))
    XYt = xcov(Xs, Y)
    m = size(XYt, 2)
#     println(size(XtY))
    res = zeros(m)
    nrms2 = mean(abs2(Y), 2)
    for kk=1:m
        val = (XYt[:,kk]'*(XcorrFact\XYt[:,kk]))[1]/nrms2[kk]
        @assert imag(val) < 1e-10
        res[kk] = real(val)
    end
    res, XcorrFact
end

function capacity(Y::Vector, Xs, XcorrFact=nothing)
    res, XcF = capacity(reshape(Y, 1, size(Y, 1)), Xs, XcorrFact)
    res[1], XcF
end

function getYs(data, maxdeg, depth)
    djks = enumerate_polys(depth, maxdeg)
    Y = legendre_prods(data, djks)
    Y, djks
end

function eval_capacities(os, Y, chunkfact, spd; discard_imag=false, Xmatrix=true)
    if discard_imag
        xdim = size(os,2)*spd
        x = real(os.')[:]
    else
        xdim = 2*size(os,2)*spd
        x = reinterpret(Float64, (os.')[:])
    end
    X = chunkwise(x, xdim*chunkfact, xdim)
    if Xmatrix
        X = Matrix(X)
    end
    capacities,_ = capacity(Y, X)
    capacities
end

function find_threshold(caps; maxiter=20, rtol=1e-2)
    converged = false
    lcaps = log(caps)
    iter = 0
    th0 = 1
    lth = log(th0)
    while !converged && iter < maxiter
        lmselected = lcaps[lcaps .<= lth]
        lm = mean(lmselected)
        lstd = stdm(lmselected, lm)
        lthdiff = (lm+3lstd-lth)*.8
        lth += lthdiff
        converged = abs(lthdiff/lth) < rtol
        iter += 1
    end
    println("Converged after $iter iterations")
    exp(lth)
end

immutable CapAnalysis
  caps
  djks
  cumcap
  lincap
  bydegmemory
  thresh
  degs
  mems
  udegs
  umems
end
  


function analyze_capacities(caps, djks, extend; plot_analysis=true)
  thresh = find_threshold(caps)
  degs = sum(djks, 1)[:]
  ccap = sum(caps[caps .> thresh])
  lccap = sum(caps[(caps .> thresh) & (degs .== 1)])
  println("capacity: ", ccap)
  println("linear capacity: ", lccap)
  println("capacity (normalized): ", ccap/extend)
  println("linear capacity (normalized): ", lccap/extend)
  
  spans = zeros(degs)
  for k=1:size(djks,2)
      non_zeros = find(djks[:,k] .> 0)
      spans[k] = non_zeros[end]-non_zeros[1]+1
  end
  uspans = unique(spans)
  udegs = unique(degs)
  capsbydegspan = zeros(length(udegs), length(uspans))
  for (kk, span)=enumerate(uspans)
      spsel = (spans.==span)&(caps.>thresh)
      for (jj, deg)=enumerate(udegs)
          degsel = (degs.==deg)
          capsbydegspan[jj,kk] = sum(caps[spsel&degsel])
      end
  end
  
  ca = CapAnalysis(caps, djks, ccap, lccap, capsbydegspan, thresh, degs, spans, udegs, uspans)  
  if plot_analysis
    plot_caps_by_j(ca)
    plot_caps_hist(ca)
    plot_caps_by_mem_deg(ca)
  end
  ca
end

function plot_caps_by_j(ca::CapAnalysis)
  f = figure(1)
  subplot(211)
  title("Information capacity per basis function", size=18)
  hlines([ca.thresh], xlim()..., linestyle="dashed")
  semilogy(ca.caps, "o")
  ylabel(L"Capacity $C[X,y_j]$", size=18)
  subplot(212)
  ylabel(L"Degree $d[y_j]$", size=18)
  plot(ca.degs)
  ylim(0, maximum(ca.degs)+1)
  xlabel(L"Function index $j$", size=18)
  td = tempdir()
  fname = "$td/caps_by_j.pdf"
  savefig("$td/caps_by_j.pdf")
  println("Saved figure at $fname")

  f
end  

function plot_caps_hist(ca::CapAnalysis)
  f = figure(2)
  plt[:hist](log(ca.caps), bins=51)
  vlines([log(ca.thresh)], ylim()..., linestyle="dashed", label="Threshold")
  xlabel(L"Capacity $C[X,y_j]$", size=18)
  ylabel("Counts", size=18)
  title("Histogram of basis capacities", size=18)
  td = tempdir()
  fname = "$td/caps_hist.pdf"
  savefig("$td/caps_hist.pdf")
  println("Saved figure at $fname")
  f
end

function plot_caps_by_mem_deg(ca::CapAnalysis)
  f = figure(3, figsize=(8,6))
  imshow(ca.bydegmemory, cmap="copper", interpolation="nearest")
  xt = 0:maximum(ca.umems)-1
  yt = 0:maximum(ca.udegs)-1
  xticks(xt, xt+1)
  yticks(yt, yt+1)
  xlabel(L"Memory $\Delta t[y_j]$", size=18)
  ylabel(L"Degree $d[y_j]$", size=18)
  title("Capacity vs. memory and degree", size=18)
  grid(false)
  colorbar()
  td = tempdir()
  fname = "$td/caps_by_mem_deg.pdf"
  savefig("$td/caps_by_mem_deg.pdf")
  println("Saved figure at $fname")

  f
end

function sample(N; seed=0)
    if seed !=0
        srand(seed)
    end
    2rand(N)-1
end



function bisect_int(fn, x0::Int, x1::Int; rtol=1e-1)
    x1 >= x0 || error("Malformed range")

    fx0 = fn(x0)
    fx1 = (x1 > x0) ? fn(x1) : fx0
    
    while true    

        xmid = floor(Int, (x0+x1)/2)
        fxmid = (xmid > x0) ? fn(xmid) : fx0
        
        if fxmid < min(fx0, fx1) && abs(1-fxmid/min(fx0, fx1)) > rtol
#             println("Non-concave, recursing.")
            x0sol, fx0sol = bisect_int(fn, x0, xmid)
            x1sol, fx1sol = bisect_int(fn, xmid, x1)
            if fx0sol < fx1sol
                return x1sol, fx1sol
            else
                return x0sol, fx0sol
            end
        end
#         println("$x0 --> $fx0")
        println("$xmid --> $fxmid")
#         println("$x1 --> $fx1")

        if fx0 < fx1
            x0 = xmid
            fx0 = fxmid
        else
            x1 = xmid
            fx1 = fxmid
        end
        
        if x1-x0 <= 1
            break
        end
    end
    if fx0 > fx1
        return x0, fx0
    else
        return x1, fx1
    end   
end
    
    

function optimize_offset(outputs, Ys, extend, spd, kminmax=(nothing, nothing))
    if kminmax == (nothing, nothing)
        kminmax = (1, 4spd)
    end
    kmin, kmax = kminmax
        
    f = kk -> sum(eval_capacities(sub(outputs, kk:size(outputs,1), 1:size(outputs,2)), Ys, extend, spd))
    kkmax, capmax = bisect_int(f, kmin, kmax)
end
        
    
