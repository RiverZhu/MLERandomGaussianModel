% bisection method

function [x_result,n] = bisect(f,a,b,delta)
n = 0;
if f(b) == 0
    x_result = b;
end
while f(a) * f(b) < 0 || isnan(f(a) * f(b))

    x_0 = (a+b)/2; 

    if f(x_0) == 0 
        x_result = x_0;
        break;
    end
    if f(a) * f(x_0) < 0 || isnan(f(a) * f(b))
        b = x_0;
    else 
        a = x_0;
    end
    n = n+1;
    
    if abs(a - b) < delta 
        x_result = x_0;
        n = n+1;
        break;
    end
    
    if (a+b)/2-x_0 == 0 
        x_result = x_0;
        break;
    end
end


