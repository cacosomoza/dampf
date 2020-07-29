function A_normalized = MPO_canonization(A,form,N,Nbe2)

% Canonical normalization of the A-matrices of the MPO: A = {A^(1),A^(2),...,A^(N)}

d = size(A{1},3); % Local physical dimension (trunctation number Nb of oscillators)

switch form
    
    % Left canonization: U. Schollwock, Annals of Physics 326,96-192 (2011). Page 129: 4.4.1 Generation of a left canonical MPS
    % sum_{sigma_l} Adag^{sigma_l} A^{sigma_l} = Id
    case 'left' 
        
        % Go through every N matrix. The A{i} matrix becomes normalized with r1 rows 
        for i1 = 1:N-1 
            [A{i1}, A_remainder_temp] = sweep_canonization( A{i1}, form, size(A{i1},1), size(A{i1},2), size(A{i1},3)) ;

            for i2 = 1:Nbe2(i1+1)                        
                A{i1+1}(:,:,i2) = A_remainder_temp*A{i1+1}(:,:,i2);                                                
            end
            
        end
        
    case 'right'
        
        for i1 = N:-1:2    
            [A{i1}, A_remainder_temp] = sweep_canonization( A{i1},form,size(A{i1},1),size(A{i1},2),size(A{i1},3) );
            for i2 = 1:d
                A{i1-1}(:,:,i2) = A{i1-1}(:,:,i2)*A_remainder_temp;
            end
        end
        
end

A_normalized = A;