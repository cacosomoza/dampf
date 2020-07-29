function [A_out, truncerr, chi_vec] = MPO_compression_adap(A_in, comp_tol, Nbe2, chimax)
    
    % Left canonization: U. Schollwock, Annals of Physics 326,96-192 (2011). Page 129: 4.4.1 Generation of a left canonical MPS
    % Right truncation:  U. Schollwock, Annals of Physics 326,96-192 (2011). Page 130: 4.5.1 Compressing a matrix product state by SVD

    % Returns number of A-matrices. A_in = A{j,k} which consists of N arrays of dim (1,chi) X (1,chi) X M^2
    M = length(A_in); 
    % Bring A-matrices into left-canonical form (without truncation)
    A_out = MPO_canonization(A_in, 'left', M, Nbe2);
                
    % Truncate A-matrices
    [A_out, truncerr, chi_vec] = MPO_truncation_adap(A_out, 'right', M, comp_tol, Nbe2, chimax);
        
end
    