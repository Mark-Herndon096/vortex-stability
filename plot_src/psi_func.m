function [return_val] = psi_func(beta)
	return_val = (beta^2)*besselk(0,beta) + beta*besselk(1,beta);
