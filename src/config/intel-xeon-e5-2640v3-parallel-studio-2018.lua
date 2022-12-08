--
-- Outer block size for the two stage schemes. These values are used by the
-- _2STAGE routines.
--

-- TGCSYLV in double precision
function tgcsylv_double_2stage(m,n)
    return 2048
end

-- TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_2stage(m,n)
    return 2048
end

-- TGCSYLV in single precision
function tgcsylv_single_2stage(m,n)
    return 2048
end

-- TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_2stage(m,n)
    return 2048
end

-- TGSYLV in double precision
function tgsylv_double_2stage(m,n)
    return 2048
end

-- TGSYLV in single precision
function tgsylv_single_2stage(m,n)
    return 2048
end

-- TRSYLV in double precision
function trsylv_double_2stage(m,n)
    return 2048
end

-- TRSYLV in single precision
function trsylv_single_2stage(m,n)
    return 2048
end

-- TRSYLV2 in double precision
function trsylv2_double_2stage(m,n)
    return 2048
end

-- TRSYLV2 in single precision
function trsylv2_single_2stage(m,n)
    return 2048
end

-- TGLYAP for double precision
function tglyap_double_2stage(m)
    return 2048
end

-- TGLYAP for single precision precision
function tglyap_single_2stage(m)
    return 2048
end

-- TRLYAP for double precision
function trlyap_double_2stage(m)
    return 2048
end

-- TRLYAP for single precision precision
function trlyap_single_2stage(m)
    return 2048
end

-- TGSTEIN for double precision
function tgstein_double_2stage(m)
    return 2048
end

-- TGSTEIN for single precision precision
function tgstein_single_2stage(m)
    return 2048
end

--
-- Block-sizes for the Level-3 and DAG schemes. These values are used by
-- routines with the names '_L3' or '_DAG'. The function with the suffix '_mb'
-- returns the block-size with respect to the number of rows of the solution.
-- The function suffixed with '_nb' returns the block-size with respect to the
-- number of columns of the solution. In case of symmetric equations, like
-- the Lypanov or the Stein equation. only the function with the suffix '_mb'
-- is used.

-- MB for TGCSYLV in double precision
function tgcsylv_double_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- NB for TGCSYLV in double precision
function tgcsylv_double_nb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- NB for TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_nb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TGCSYLV in single precision
function tgcsylv_single_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- NB for TGCSYLV in single precision
function tgcsylv_single_nb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- NB for TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_nb(m,n)


end

-- MB for TGSYLV in double precision
function tgsylv_double_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    else
        return 96
    end
end


-- NB for TGSYLV in double precision
function tgsylv_double_nb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    elseif ( math.min(m.n) <= 2816 ) then
        return 96
    else
        return 128
    end
end

-- MB for TGSYLV in single precision
function tgsylv_single_mb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    elseif ( math.min(m.n) <= 2816 ) then
        return 96
    else
        return 128
    end
end

-- NB for TGSYLV in single precision
function tgsylv_single_nb(m,n)
    if ( math.min(m,n) <= 2048 ) then
        return 64
    elseif ( math.min(m.n) <= 2816 ) then
        return 96
    else
        return 128
    end
end

-- MB for TRSYLV in double precision
function trsylv_double_mb(m,n)
    if ( math.min(m,n) <= 512 ) then
        return 32
    elseif ( math.min(m.n) <= 2304 ) then
        return 64
    else
        return 96
    end
end


-- NB for TRSYLV in double precision
function trsylv_double_nb(m,n)
    if ( math.min(m,n) <= 512 ) then
        return 32
    elseif ( math.min(m.n) <= 2304 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRSYLV in single precision
function trsylv_single_mb(m,n)
    if ( math.min(m,n) <= 512 ) then
        return 32
    elseif ( math.min(m.n) <= 2304 ) then
        return 64
    else
        return 96
    end
end

-- NB for TRSYLV in single precision
function trsylv_single_nb(m,n)
    if ( math.min(m,n) <= 512 ) then
        return 32
    elseif ( math.min(m.n) <= 2304 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRSYLV2 in double precision
function trsylv2_double_mb(m,n)
    if ( math.min(m,n) <= 1792 ) then
        return 64
    elseif ( math.min(m.n) <= 3072 ) then
        return 96
    else
        return 128
    end
end

-- NB for TRSYLV2 in double precision
function trsylv2_double_nb(m,n)
    if ( math.min(m,n) <= 1792 ) then
        return 64
    elseif ( math.min(m.n) <= 3072 ) then
        return 96
    else
        return 128
    end
end

-- MB for TRSYLV2 in single precision
function trsylv2_single_mb(m,n)
    if ( math.min(m,n) <= 1792 ) then
        return 64
    elseif ( math.min(m.n) <= 3072 ) then
        return 96
    else
        return 128
    end
end

-- NB for TRSYLV2 in single precision
function trsylv2_single_nb(m,n)
    if ( math.min(m,n) <= 1792 ) then
        return 64
    elseif ( math.min(m.n) <= 3072 ) then
        return 96
    else
        return 128
    end
end

-- MB for TGLYAP in double precision
function tglyap_double_mb(m)
    if ( m <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TGLYAP in single precision
function tglyap_single_mb(m)
    if ( m <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRLYAP in double precision
function trlyap_double_mb(m)
    if ( m <= 3072 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRLYAP in single precision
function trlyap_single_mb(m)
    if ( m <= 3072) then
        return 64
    else
        return 96
    end
end

-- MB for TGSTEIN in double precision
function tgstein_double_mb(m)
    if ( m <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TGSTEIN in single precision
function tgstein_single_mb(m)
    if ( m <= 2048 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRSTEIN in double precision
function trstein_double_mb(m)
    if ( m <= 3072 ) then
        return 64
    else
        return 96
    end
end

-- MB for TRSTEIN in single precision
function trstein_single_mb(m)
    if ( m <= 3072 ) then
        return 64
    else
        return 96
    end
end


