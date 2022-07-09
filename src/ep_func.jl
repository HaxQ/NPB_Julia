
function randlc(x, a)
    i246m1 = 0x00003FFFFFFFFFFF
    d2m46 = 0.5^46

    Lx = trunc(Int64, x)
    La = trunc(Int64, a)

    Lx = (Lx * La) & i246m1
    x = Float64(Lx)
    return x, (d2m46 * Float64(Lx))
end

function vranlc(n, x, a, y)
    i246m1 = 0x00003FFFFFFFFFFF
    d2m46 = 0.5^46

    Lx = trunc(Int64, x)
    La = trunc(Int64, a)

    for i = 1:n
        Lx = (Lx * La) & i246m1
        y[i] = d2m46 * Float64(Lx)
    end
    x = Float64(Lx)
    return x
end

function verify(m, sx, sy, gc)
    epsilon = 1.0E-08
    verified = true
    if m == 28
        class = 'A'
        sx_verify_value = 1.682235632304711E+08
        sy_verify_value = 1.682195123368299E+08
        gc_verify_value = 210832767.0
    elseif m == 30
        class = 'B'
        sx_verify_value = 6.728927543423024E+08
        sy_verify_value = 6.728951822504275E+08
        gc_verify_value = 843345606.0
    elseif m == 32
        class = 'C'
        sx_verify_value = 2.691444083862931E+09
        sy_verify_value = 2.691519118724585E+09
        gc_verify_value = 3373275903.0
    end

    if (verified)
        sx_err = abs((sx - sx_verify_value)/sx_verify_value)
        sy_err = abs((sy - sy_verify_value)/sy_verify_value)
        if (isnan(sx_err) || isnan(sy_err))
        verified = false
        else
        verified = ((sx_err <= epsilon) && (sy_err <= epsilon))
        end
    end

    if (verified)
        gc_err = abs((gc - gc_verify_value)/gc_verify_value)
        if (isnan(gc_err) || (gc_err > epsilon))
        verified = false
        end
    end
    return verified
end

function compute_gaussian_deviates(nk, x, qq, sx, sy)
    #for i = 1:nk 
    @fastmath @inbounds @simd for i = 1:nk
        x1 = 2.0 * x[2*i - 1] - 1.0
        x2 = 2.0 * x[2*i] - 1.0
        t1 = x1^2 + x2^2
        if t1 <= 1.0
            t2 = sqrt(-2.0 * log(t1) / t1)
            t3 = abs(x1 * t2)
            t4 = abs(x2 * t2)
            l = trunc(UInt8, max(t3, t4))
            qq[l+1] += 1.0
            sx += t3
            sy += t4
        end
    end
    
    return sx, sy
end

function generate_random_numbers(kk, nk, t1, t2, a, x)
    for i = 1:100
        ik = div(kk, 2)
        if (2 * ik) != kk
            t1, t3 = randlc(t1, t2)
        end
        if ik == 0
            break
        end
        t2, t3 = randlc(t2, t2)
        kk = ik
    end
    t1 = vranlc(2 * nk, t1, a, x)
end