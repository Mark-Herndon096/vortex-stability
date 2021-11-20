function [return_val] = chi_func(beta)
	return_val = beta*besselk(1,beta);
