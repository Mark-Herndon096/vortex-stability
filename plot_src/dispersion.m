function [return_val] = dispersion(beta)
    global ka_val;
	J0 = besselj(0,beta);
	J1 = besselj(1,beta);
	J2 = besselj(2,beta);
	K0 = besselk(0,ka_val);
	K1 = besselk(1,ka_val);
	K2 = besselk(2,ka_val);

	J1_p =  (J0 - J2)/2.0;
	K1_p = -(K0 + K2)/2.0;
	
	return_val = (1.0/beta)*(J1_p/J1) + K1_p/(ka_val*K1) +  ...
				 sqrt(beta^2 + (ka_val)^2)/(ka_val*beta^2);
