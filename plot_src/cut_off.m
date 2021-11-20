function [omega] = cut_off(delta)
    delta = 0.642*delta; % For a rankine vortex
	t1 = (cos(delta) - 1.0)/delta^2;
	t2 =  sin(delta)/delta;
	t3 = cosint(delta);
	
	omega = (1/2)*(t1 + t2 - t3);
