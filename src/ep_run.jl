using Printf

function ep(class)
    # data
    if class == 'A'
        m = 28
    elseif class == 'B'
        m = 30
    elseif class == 'C'
        m = 32
    end

    mk = 16
    mm = m - mk
    nn = 2^mm
    nk = 2^mk
    nq = 10

    a = 1220703125.0
    s = 271828183.0

    ## storage
    x = Array{Float64, 1}(undef, 2*nk)

    #since lower bound is not variable in Julia, may have to tweak index math
    qq = Array{Float64, 1}(undef, nq)
    #qq = MVector{nq,Float64}
    q = Array{Float64, 1}(undef, nq)
    #q = MVector{nq,Float64}

    # end data

    #initialization

    println("\n\n NAS Parallel Benchmarks (NPB3.4-Julia) - EP Benchmark\n")
    size = 2.0^(m+1)
    println(" Number of random numbers generated: ", Int64(size))

    np = nn

    Mops = log(sqrt(abs(max(1.0,1.0))))

    for i = 1:2*nk
        x[i] = -1.0e99
    end

    # start timer

    t1 = a
    t1 = vranlc(0, t1, a, x)
    # Compute AN = A ^ (2 * NK) (mod 2^46)

    t1 = a
    t2 = 0.0
    for i = 1:mk+1
        t1, t2 = randlc(t1, t1)
    end

    an = t1
    tt = s
    gc = 0.0
    sx = 0.0
    sy = 0.0
    for i = 1:nq
        q[i] = 0.0
    end

    k_offset = -1

    for i = 1:nq
        qq[i] = 0.0
    end
    
    timerand = 0.0
    timecompute = 0.0
    for k = 1:np
        kk = k_offset + k
        t1 = s
        t2 = an

        timerand += @elapsed generate_random_numbers(kk, nk, s, an, a, x)
        timecompute += @elapsed sx, sy = compute_gaussian_deviates(nk, x, qq, sx, sy)
    end

    for i = 1:nq
        q[i] += qq[i]
    end

    for i = 1:nq
        gc += q[i]
    end

    # stop timer

    verified = verify(m, sx, sy, gc)
    println("sx: $sx sy: $sy gc: $gc")
    if(verified)
        @printf("Verification = SUCCESSFUL\n")
    else
        @printf("Verification = FAILED\n")
    end

    Mops = 2.0^(m+1)/(timerand + timecompute)/1000000.0
    
    @printf("Mops/s total: %0.2f\n", Mops)
    println()
   
    @printf("Gaussian pairs (seconds): %0.3f\n", timecompute)
    @printf("Random numbers (seconds): %0.3f\n", timerand)
end

#totaltime = @elapsed main_prog()
#@printf("Total time: %0.3f\n", totaltime)