function [root] = root_find(my_func, a, b, tol)
    x1 = a; x2 = b;
    
    for i = 1:60
        mid  = (x1 + x2)/2.0;
        product = my_func(x1)*my_func(mid);
        if  product <= 0.0
            x2  = mid;
            mid = (x1 + x2)/2.0;
            dx  = abs(x2 - x1);
        elseif  product > 0.0 
            x1  = mid;
            mid = (x1 + x2)/2.0;
            dx  = abs(x1 - x2);
        end
        
        if  dx <= tol
            break;
        end
       
    end
    
    root = mid;
   % fprintf('Root approximated at %16.16f\n',root)
