function [A_normalized, A_remainder] = sweep_canonization(A,form,D_L,D_R,d)

% transform MP parametrization to canonical parametrization on a single site.
% A = A^(n), with site-index n.
%---------------------------------------------------------
% Example: 'left'
% transform to left-canonical form \sum_s A_s' * A_s = id
% (1) A_stored_vertical = [A_1;A_2;...;A_d]
% (2) compute SVD: [U, S, V] = svd(A_stored_vertical,'econ')
% (3) A_left_normalized = U, with remainder S * V' (A_s = A_l_n_s * S * V')
%
% Example: 'right'
% similar :)

switch form
    
    case 'left'
                
        if D_L*d >= D_R 
                                                
            A_stored_vertical = reshape(permute(A,[1,3,2]),D_L*d,D_R);                        
            [U,R] = qr(A_stored_vertical,0);                    
            A_normalized = permute(reshape(U,D_L,d,D_R),[1,3,2]);
            A_remainder = R;
                        
        else            
            
            U = zeros(D_L*d,D_R); 
            R = zeros(D_R,D_R);   
            A_stored_vertical = reshape(permute(A,[1,3,2]),D_L*d,D_R);                                
            [U(1:D_L*d,1:D_L*d), R(1:D_L*d,1:D_R)] = qr(A_stored_vertical,0);                                                
            A_normalized = permute( reshape(U, D_L, d, D_R), [1,3,2] );                                               
            A_remainder = R;        
        end
   
    case 'right'
        
        if D_L <= D_R*d
                    
            A_stored_horizontal = reshape(A,D_L,D_R*d);
            [U, S, V] = svd(A_stored_horizontal,'econ');
            A_normalized = reshape(V',D_L,D_R,d);
            A_remainder = U*S;
        
        else

            U = zeros(D_L,D_L);
            S = zeros(D_L,D_L);
            V =zeros(D_R*d,D_L);
            
            A_stored_horizontal = reshape(A,D_L,D_R*d);
            [U(1:D_L,1:D_R*d), S(1:D_R*d,1:D_R*d), V(1:D_R*d,1:D_R*d)] = svd(A_stored_horizontal,'econ');
            A_normalized = reshape(V',D_L,D_R,d);
            A_remainder = U*S;
            
        end
        
    otherwise
        
        warning('Unexpected normalization type.')
  
end