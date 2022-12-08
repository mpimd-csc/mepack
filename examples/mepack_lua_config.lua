--
-- Example MEPACK config file
-- Copyright (C) Martin Koehler, 2017-2022
--

--
-- Blocksizes for double precision routines
--

-- TRSYLV row block size (MB)
function trsylv_double_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV column block size (NB)
function trsylv_double_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV2 row block size (MB)
function trsylv2_double_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV2 column block size (NB)
function trsylv2_double_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSYLV row block size (MB)
function tgsylv_double_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSYLV column block size (NB)
function tgsylv_double_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV row block size (MB)
function tgcsylv_double_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV column block size (NB)
function tgcsylv_double_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV_DUAL row block size (MB)
function tgcsylv_dual_double_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV_DUAL column block size (NB)
function tgcsylv_dual_double_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end


-- TRLYAP block size (MB)
function trlyap_double_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSTEIN block size (MB)
function trstein_double_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGLYAP block size (MB)
function tglyap_double_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSTEIN block size (MB)
function tgstein_double_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

--
-- Two stage block sizes
--

-- TRLYAP two stage block size
function trlyap_double_2stage(m)
    return 3072
end

-- TGLYAP two stage block size
function tglyap_double_2stage(m)
    return 3072
end

-- TRSTEIN two stage block size
function trstein_double_2stage(m)
    return 3072
end

-- TGSTEIN two stage block size
function tgstein_double_2stage(m)
    return 3072
end

-- TRSYLV two stage block size
function trsylv_double_2stage(m,n)
    return 3072
end

-- TRSYLV two stage block size
function trsylv2_double_2stage(m,n)
    return 3072
end

-- TGSYLV two stage block size
function tgsylv_double_2stage(m,n)
    return 3072
end

-- TGCSYLV two stage block size
function tgcsylv_double_2stage(m,n)
    return 3072
end

-- TGCSYLV_DUAL two stage block size
function tgcsylv_dual_double_2stage(m,n)
    return 3072
end










--
-- Blocksizes for single precision routines
--

-- TRSYLV row block size (MB)
function trsylv_single_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV column block size (NB)
function trsylv_single_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV2 row block size (MB)
function trsylv2_single_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSYLV2 column block size (NB)
function trsylv2_single_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSYLV row block size (MB)
function tgsylv_single_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSYLV column block size (NB)
function tgsylv_single_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV row block size (MB)
function tgcsylv_single_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV column block size (NB)
function tgcsylv_single_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV_DUAL row block size (MB)
function tgcsylv_dual_single_mb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGCSYLV_DUAL column block size (NB)
function tgcsylv_dual_single_nb(m,n)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end


-- TRLYAP block size (MB)
function trlyap_single_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TRSTEIN block size (MB)
function trstein_single_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGLYAP block size (MB)
function tglyap_single_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

-- TGSTEIN block size (MB)
function tgstein_single_mb(m)
    if ( m <= 512 and n <= 512 ) then
        return 32
    else
        return 64
    end
end

--
-- Two stage block sizes
--

-- TRLYAP two stage block size
function trlyap_single_2stage(m)
    return 3072
end

-- TGLYAP two stage block size
function tglyap_single_2stage(m)
    return 3072
end

-- TRSTEIN two stage block size
function trstein_single_2stage(m)
    return 3072
end

-- TGSTEIN two stage block size
function tgstein_single_2stage(m)
    return 3072
end

-- TRSYLV two stage block size
function trsylv_single_2stage(m,n)
    return 3072
end

-- TRSYLV two stage block size
function trsylv2_single_2stage(m,n)
    return 3072
end

-- TGSYLV two stage block size
function tgsylv_single_2stage(m,n)
    return 3072
end

-- TGCSYLV two stage block size
function tgcsylv_single_2stage(m,n)
    return 3072
end









