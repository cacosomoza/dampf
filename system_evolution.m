function [A, truncerr, chivals] = system_evolution( A, it, Aid, epsilon, N, Q, Us, oide, Nbe, sgntrnsp, comp_tol, chivals, parallel, chimax )
   

    % Blank MPO
    BlankMPO = cell(1,N);
    for k = 1:N*Q 
        BlankMPO{k} = zeros(1,1,Nbe(k)^2);
    end

    A_temp = cell(numel(Aid),1);  
    truncerr = 0; 
    
    if parallel
                                
            parfor i = 1:numel(Aid)

                [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables                       
                A_temp{i} = BlankMPO;          % Auxiliary blank MPO (with 1 X 1 matrices)            
                trer = 0;   
                trertemp=0;

                for m = 1:numel(Aid) % for every Omn        

                    [s,t] = ind2sub([N,N],Aid(m)); % extract the m,n indeces
                    %chi_vec = ones(1,N*Q-1);       % bond results of a given sum term
                    
                    if abs( Us(j, k, s, t) ) > epsilon

                        A_temp{i} = MPO_sum( A_temp{i}, MPO_scalar( Us(j, k, s, t), A{m}) , Nbe.^2);
                        % The (n,m) term needs is extracted from conjugate A{m,n}
                        if s ~= t
                            A_temp{i} = MPO_sum( A_temp{i}, MPO_scalar( Us(j, k, t, s), MPO_conjugate( A{m}, sgntrnsp, Nbe.^2, oide )), Nbe.^2);
                        end
                        [ A_temp{i}, trertemp, chi_vec ] = MPO_compression_adap(A_temp{i}, comp_tol , Nbe.^2, chimax);                                                                
                    end
                    
                    
                    trer = trer + trertemp;                    
                end                            
                if j ~= k
                    trer = 2*trer;
                end                
                
                
                chivals{i}(it,:) = chi_vec ;
                truncerr = truncerr + trer;                                                                     
            end
            A = A_temp;    
    else
                             
            for i = 1:numel(Aid)

                [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables                       
                A_temp{i} = BlankMPO;          % Auxiliary blank MPO (with 1 X 1 matrices)            
                trer = 0;   
                trertemp=0;

                for m = 1:numel(Aid) % for every Omn        

                    [s,t] = ind2sub([N,N],Aid(m)); % extract the m,n indeces
                    chi_vec = ones(1,N*Q-1);

                    if abs(Us(j, k, s, t)) > epsilon

                        A_temp{i} = MPO_sum( A_temp{i}, MPO_scalar( Us(j, k, s, t), A{m}) , Nbe.^2);
                        % The (n,m) term needs is extracted from conjugate A{m,n}
                        if s ~= t
                            A_temp{i} = MPO_sum( A_temp{i}, MPO_scalar( Us(j, k, t, s), MPO_conjugate( A{m}, sgntrnsp, Nbe.^2, oide )), Nbe.^2);
                        end
                        [A_temp{i}, trertemp, chi_vec] = MPO_compression_adap(A_temp{i}, comp_tol , Nbe.^2, chimax);                                                                

                    end                                        
                    trer = trer + trertemp;                   
                end                            
                if j ~= k
                    trer = 2*trer;
                end                
                
                chivals{i}(it,:) = chi_vec ;
                truncerr = truncerr + trer;                                                                         
            end
            A = A_temp;            
    end
   

end