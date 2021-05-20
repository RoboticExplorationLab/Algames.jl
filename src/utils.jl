################################################################################
# Helpers
################################################################################

function add2sub(v::SubArray, e)
	v .+= e
	return nothing
end

function add2sub(v::SubArray, e::Diagonal{T,SVector{n,T}}) where {n,T}
	for i = 1:n
		v[i,i] += e[i,i]
	end
	return nothing
end

function addI2sub(v::SubArray, e)
	n = size(v)[1]
	for i = 1:n
		v[i,i] += e
	end
	return nothing
end

function sparse_zero!(spm::SparseMatrixCSC)
	n = length(spm.nzval)
	for i = 1:n
		spm.nzval[i] = 0.0
	end
	return nothing
end

################################################################################
# Printers
################################################################################

function display_solver_header()
	@printf(
		"%-3s %-2s %-2s %-6s %-6s %-6s \n",
		"out",
		"in",
		"α",
		"Δ",
		"res",
		"reg",
		)
	return nothing
end

function display_solver_data(k, l, j, Δ, res_norm, reg)#, condi, loss, val_scale, jac_scale)
	@printf(
		"%-3s %-2s %-2s %-6s %-6s %-6s \n",
		k,
		l,
		j,
		@sprintf("%.0e", Δ),
		@sprintf("%.0e", res_norm),
		@sprintf("%.0e", reg.x),
		)
	return nothing
end

function scn(a::Number; digits::Int=1)
	@assert digits >= 0
    # a = m x 10^e
    if a == 0
        e = 0
        m = 0.0
    else
        e = Int(floor(log(abs(a))/log(10)))
        m = a*exp(-e*log(10))
    end
    m = round(m, digits=digits)
    if digits == 0
        m = Int(floor(m))
		strm = string(m)
	else
		strm = string(m)
		is_neg = m < 0.
		strm = strm*"0"^abs(2+digits+is_neg-length(strm))
    end
    sgn = a >= 0 ? " " : ""
    sgne = e >= 0 ? "+" : ""
    return "$sgn$(strm)e$sgne$e"
end
