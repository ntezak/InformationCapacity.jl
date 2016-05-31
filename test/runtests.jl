using InformationCapacity
using Base.Test

using Jacobi

function test_chunked(;
    Ns=[100,232,501], 
    cls=[25,30],
    css=[10,25]
    )
    for N=Ns
        for cl=cls
            for cs=css
                x = randn(N)
                X = x2X(x, cl, cs)
                xc = chunkwise(x, cl, cs)
                xcM = Matrix(xc)
                @test norm(xcM-X) < 1e-5
            end
        end
    end
end

function test_xcov(;
    Ns=[100,232,501], 
    cls=[25,30],
    css=[10,25])
    for N=Ns
        for cl=cls
            for cs=css
                x = randn(N)
                y = randn(2N)
                xc = chunkwise(x, cl, cs)
                yc = chunkwise(y, 2cl, 2cs)                
                xy = xcov(xc, yc)
                xy2 = (Matrix(xc)*Matrix(yc)')/size(xc,2)
                @test norm(xy-xy2) < 1e-5
            end
        end
    end
end

function test_enumerate_polys()
    djks1 = enumerate_polys_bydeg(2,1)
    djks2 = enumerate_polys_bydeg(2,2)
    @test norm(djks1 - [1 0; 0 1]) < 1e-5
    @test norm(djks2 - [2 1 0; 0 1 2]) < 1e-5
    djks12 = enumerate_polys(2,2)
    @test norm(djks12-[djks1 djks2]) < 1e-5
end

function test_legendre_prods(N=10)
    u = 2rand(N)-1
    djks1 = enumerate_polys_bydeg(2,1)
    djks2 = enumerate_polys_bydeg(2,2)
    y1 = legendre(u, 1)
    y2 = legendre(u, 2)
    y11 = y1[1:N-1]
    y12 = y1[2:N]
    
    y21 = y2[1:N-1]
    y22 = y2[2:N]
    
    @test norm(legendre_prods(u, djks1)-[y11'; y12']) < 1e-5
    @test norm(legendre_prods(u, djks2)-[y21'; (y11.*y12)'; y22']) < 1e-5
end

function test_capacity(N=100_000)
    u = 2*rand(N)-1
    djks12 = enumerate_polys(2,2)
    djks3 = enumerate_polys_bydeg(2,3)    
    y12 = legendre_prods(u, djks12)
    y3 = legendre_prods(u, djks3)
    x12 = y12
    c12, xacov = capacity(y12, x12)
    @test norm(c12-ones(size(c12)...)) < 1e-5
    c3, xacov = capacity(y3, x12)
    @test norm(c3) < 1e-3
    
end


function test_full_analysis(N0=5000)
  
  chi = 1.
  beta = 10.
  sigma = 50.

  Δ0 = .5
  u0s = 0
  u0s = 20
  u0 = (-1-1im)*u0s


  spd0 = 16


  u = sample(N0)
  hmax = .005
  println("Simulating DOPO reservor by modulating pump")
  @time us, os, nle = DOPO.simulate_reservoir2(u, chi, beta, Δ0, sigma, u0; 
    hmax=hmax, spd=spd0);
      
  cf0 = 4
  md0 = 8
  depth0 = 5
  extend=1
  
  println("Computing target output basis functions")
  @time Ys, djks = getYs(u, md0, depth0);
  
  println("Optimizing data offset")
  @time kkopt, copt = optimize_offset(os, Ys, extend, spd0)
  
  @test copt >= 10.
  
  oskk = sub(os, kkopt:size(os,1), 1:size(os, 2))
  
  println("Computing capacities")
  @time caps = eval_capacities(oskk, Ys, extend, spd0)
  
  println("Performing analysis")
  ca = analyze_capacities(caps, djks, extend)
  @test ca.cumcap >= 10. 
end
    
test_chunked()
test_xcov()
test_enumerate_polys() 
test_legendre_prods()
test_capacity()
test_full_analysis()
