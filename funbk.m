%Function that calculates the observation vectors. 
function b_k = funbk(x,v)
    theta_t = atan(x(2)./x(1));
    r_p = (v'*x)./norm(x,2);
    b_k = [theta_t;r_p];
end
