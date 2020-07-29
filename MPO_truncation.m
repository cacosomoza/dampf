function [A_normalized, truncerr, chi_vec] = MPO_truncation(A,form,M,comp_tol,Nbe2)

%normalize MP parametrization A = {A^(1),A^(2),...,A^(N)}

d = size(A{1},3);

switch form
    
    % case 'left'  % NOT IMPLEMENTED YET
        
    %     for i1 = 1:N-1    
    %         [A{i1}, A_remainder_temp] = sweep_truncation(A{i1},form,size(A{i1},1),size(A{i1},2),size(A{i1},3));
    %         for i2 = 1:d
    %             A{i1+1}(:,:,i2) = A_remainder_temp*A{i1+1}(:,:,i2);
    %         end
    %     end
        
    % Right truncation:  U. Schollwock, Annals of Physics 326,96-192 (2011). Page 130: 4.5.1 Compressing a matrix product state by SVD
    case 'right'
        
        chi_vec = zeros(1,M-1);
        
        truncerr = 0;
        for i1 = M:-1:2    
            [A{i1}, A_remainder_temp, truncerrtemp, chi] = sweep_truncation( A{i1},form,size(A{i1},1),size(A{i1},2),size(A{i1},3), comp_tol, chimax);
            
            % caution: size of A{i} may change with truncation!                         
            truncerr = truncerr + truncerrtemp;
            for i2 = 1:Nbe2(i1-1)                
                A_temp{i1-1}(:,:,i2) = A{i1-1}(:,:,i2)*A_remainder_temp;
            end            
            A{i1-1} = A_temp{i1-1};        
            
            chi_vec(i1-1) = chi;
            
        end
end

A_normalized = A;

