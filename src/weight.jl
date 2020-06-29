function weight(PPDE, pThresh, alpha, precision)
	((PPDE >= 0) & (PPDE <= 1)) || error("PPDE has to be between 0 and 1")
	if PPDE==1
	        weight = -alpha * 2 * log(1-pThresh)
	else
	        weight = min(2 * (-log(1 - PPDE) + log(1 - pThresh)),
				-alpha * 2 * log(1 - pThresh))
	end

	return round(weight, sigdigits = precision)
end
