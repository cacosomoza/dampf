function rho = normalize_rho(rho)

    
    timesize = size(rho,1);
    N   = sqrt(size(rho,2));

    sysID=reshape(find(ones(N,N)),[N,N]);                        % 1d idx of rho_S in matrix form 
    pop_idx=diag(sysID);  
    
    for tt=1:timesize       
       rho(tt,:) = rho(tt,:)/ sum( rho(tt,pop_idx) );        
    end

end