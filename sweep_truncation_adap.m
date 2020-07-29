    function [A_normalized, A_remainder, truncerr, chi] = sweep_truncation_adap(A,form,D_L,D_R,d,comp_tol,chimax)

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

% set options
    % opts.tol = 1e-8;
    % opts.maxit = 150;

switch form
    
    case 'left' % not employed
        
        if D_L*d >= D_R
            
            A_stored_vertical = reshape(permute(A,[1,3,2]),D_L*d,D_R);
            [U, S, V] = svd(A_stored_vertical);
            A_normalized = permute(reshape(U,D_L,d,D_R),[1,3,2]);
            A_remainder = S*V';
        
        else

            U = zeros(D_L*d,D_R);
            S = zeros(D_R,D_R);
            V =zeros(D_R,D_R);
            A_stored_vertical = reshape(permute(A,[1,3,2]),D_L*d,D_R);
            [U(1:D_L*d,1:D_L*d), S(1:D_L*d,1:D_L*d), V(1:D_R,1:D_L*d)] = svd(A_stored_vertical);
            A_normalized = permute(reshape(U,D_L,d,D_R),[1,3,2]);
            A_remainder = S*V';
        
        end
                       
    % Truncated SVD
    case 'right'
        % SVD econ means that if NA<NB or viceversa, U is truncated to the minimum dim    
        if D_L <= D_R*d

            A_stored_horizontal = reshape(A,D_L,D_R*d);
            
            try                

                [U, S, V] = svd(A_stored_horizontal,'econ');           
               %[U, S, V,~] = lmsvd(A_stored_horizontal,floor(min(size(A_stored_horizontal))/2),opts);           
                
            catch
                disp('case 1: svd( ,econ) did not converge. Size of input matrix: \n')
                size( A_stored_horizontal )
                disp('Any NaN?')
                sum(sum(isnan( A_stored_horizontal )))
                disp('Any inf')
                sum(sum(isinf( A_stored_horizontal )))                
                disp('Trying svds(,flag) \n')
                [U, S, V,flag] = svds(A_stored_horizontal,1); 
                if flag == 1
                   disp('some singular values did not converge') 
                end
            end                       
            
            % Analise decay of singular values
            Svec = diag(S);
            Snorm = sum(Svec.^2); % as a vector                                    
            comp_cond = 1;
            cc=1;
            while comp_cond && (cc <= chimax)                                  
                 uptocc = sum(Svec(1:cc).^2) / Snorm;       % should be close to 1           
                 %cc, uptocc, Snorm                 
                 if uptocc < comp_tol
                     cc = cc+1;
                 else
                    comp_cond = 0; % get out                    
                 end                 
            end 
            chi = cc; % successfull chi
            %chi = min(chi,min(size(S)));

           
            A_normalized = reshape(V(:,1:chi)',chi,D_R,d);                  
            A_remainder = U(:,1:chi)*S(1:chi,1:chi);            
            if chi+1 < size(S,1)
                truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
            else
                truncerr = 0;
            end
            
        else

            U = zeros(D_L,D_L);
            S = zeros(D_L,D_L);
            V =zeros(D_R*d,D_L); 
            A_stored_horizontal = reshape(A,D_L,D_R*d);                        
            
            try
                
                [U(1:D_L,1:D_R*d), S(1:D_R*d,1:D_R*d), V(1:D_R*d,1:D_R*d)] = svd(A_stored_horizontal,'econ');                
                %[U, S, V,] = lmsvd(A_stored_horizontal,floor(min(size(A_stored_horizontal))/2),opts);           
                
            catch
                disp('case 2: svd( ,econ) did not converg. Size of input matrix: \n')
                size( A_stored_horizontal )
                disp('Any NaN?')
                sum(sum(isnan( A_stored_horizontal )))
                disp('Any inf')
                sum(sum(isinf( A_stored_horizontal )))                
                disp('Trying svds(,flag) \n')                                                                             
                
                [U, S, V,flag] = svds(A_stored_horizontal,1); 
                                
                if flag==1
                   disp('some singular values did not converge') 
                end
            end
                      
            % Analise decay of singular values
            Svec = diag(S);
            Snorm = sum(Svec.^2); % as a vector                                    
            comp_cond = 1;
            cc=1;
            while comp_cond && (cc <= chimax)                 
                 uptocc = sum(Svec(1:cc).^2) / Snorm ;                 
                 if uptocc < comp_tol
                     cc = cc+1;
                 else
                    comp_cond = 0; % get out                    
                 end                 
            end 
            chi = cc; % successfull chi
            %chi = min(chi,min(size(S)));
                        
            A_normalized = reshape(V(:,1:chi)',chi,D_R,d);
            A_remainder = U(:,1:chi)*S(1:chi,1:chi);
            
            if chi+1 < size(S,1)
                truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
            else
                truncerr = 0;
            end

        end
        
% % Randomised SVD
%     case 'right'
%         
%         if D_L <= D_R*d
%         
%             A_stored_horizontal = reshape(A,D_L,D_R*d);
%             %[U, S, V] = svd(A_stored_horizontal,'econ');
%             [U, S, V] = rsvd(A_stored_horizontal,chi);
%             
%             chi = min(chi,min(size(S)));
%             A_normalized = reshape(V(:,1:chi)',chi,D_R,d);                  
%             A_remainder = U(:,1:chi)*S(1:chi,1:chi);            
%             if chi+1 < size(S,1)
%                 truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
%             else
%                 truncerr = 0;
%             end
%             
%         else
% 
%             U = zeros(D_L,D_L);
%             S = zeros(D_L,D_L);
%             V =zeros(D_R*d,D_L);
%             A_stored_horizontal = reshape(A,D_L,D_R*d);
%             %[U(1:D_L,1:D_R*d), S(1:D_R*d,1:D_R*d), V(1:D_R*d,1:D_R*d)] = svd(A_stored_horizontal,'econ');
%             [U(1:D_L,1:D_R*d), S(1:D_R*d,1:D_R*d), V(1:D_R*d,1:D_R*d)] = rsvd(A_stored_horizontal,chi);
%             
%             chi = min(chi,min(size(S)));
%             A_normalized = reshape(V(:,1:chi)',chi,D_R,d);
% 
%             A_remainder = U(:,1:chi)*S(1:chi,1:chi);
%             
%             if chi+1 < size(S,1)
%                 truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
%             else
%                 truncerr = 0;
%             end
% 
%         end
        
    otherwise
        
        warning('Unexpected normalization type.')
  
end


%     case 'right'
%         
%         if D_L <= D_R*d
%         
%             A_stored_horizontal = reshape(A,D_L,D_R*d);
%             [U, S, V] = rsvd(A_stored_horizontal,chi);
%             
%             chi = min(chi,min(size(S)));
%             A_normalized = reshape(V(:,1:chi)',chi,D_R,d);                  
%             A_remainder = U(:,1:chi)*S(1:chi,1:chi);            
%             if chi+1 < size(S,1)
%                 truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
%             else
%                 truncerr = 0;
%             end
%             
%         else
% 
%             U = zeros(D_L,D_L);
%             S = zeros(D_L,D_L);
%             V =zeros(D_R*d,D_L);
%             A_stored_horizontal = reshape(A,D_L,D_R*d);
%             [U(1:D_L,1:D_R*d), S(1:D_R*d,1:D_R*d), V(1:D_R*d,1:D_R*d)] = svd(A_stored_horizontal,'econ');%             
%             
%             chi = min(chi,min(size(S)));
%             A_normalized = reshape(V(:,1:chi)',chi,D_R,d);
% 
%             A_remainder = U(:,1:chi)*S(1:chi,1:chi);
%             
%             if chi+1 < size(S,1)
%                 truncerr = sum(diag(S(chi+1:end,chi+1:end)).^2);
%             else
%                 truncerr = 0;
%             end
% 
%         end
