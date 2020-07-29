 
% A_in:     input MPO 
% isite:    on which site the A-matrix changes (A_i^(isite))_jk
% y:        the local linear transformation given as a matrix. Size(M^2 X M^2)   

function A_out = MPO_local_update(A_in, isite, y)
    
    % A_in{isite}: 1 X chi X M^2, chi X chi X M^2 or chi X 1 X M^2
    D1 = size(A_in{isite}, 1); % 1 or chi
    D2 = size(A_in{isite}, 2); % 1 or chi
    d  = size(A_in{isite}, 3); % Nb^2

    A_out = A_in;
            
    A_out{isite} = reshape( reshape(A_in{isite}, [], d) * y.', D1, D2, d );

end