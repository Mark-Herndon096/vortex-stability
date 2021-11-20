function [return_val] = self_induct(ka)
    global ka_val
    
    ka_val = ka;
	a = 0.01; b = 5.0; tol = 1e-14; eps = 1e-10;
	root_val = root_find(@special, a, b, tol);
	b = root_val - eps;
	root_val = root_find(@dispersion, a, b, tol);
	return_val = ((2.0*ka/sqrt(ka^2 + root_val^2)) - 1.0);
