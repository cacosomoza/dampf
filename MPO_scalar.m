function A_out = MPO_scalar(x, A_in)

    A_out = A_in;    
    A_out{1} = x*A_in{1};

end