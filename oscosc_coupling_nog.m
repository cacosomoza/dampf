function [A, bonderr] = oscosc_coupling( A, N, Q, Aid, U_osc_int, Nb, Nbe, comp_tol, parallel)

    if parallel

        bonderr = 0;
              
            parfor i=1:numel(Aid)             
                for nn = 1:N
                for bb = 1:Q-1    
                    % Ex: Q=3, M=6, isite= (1)-2, (2)-3, (4)-5, (5)-6
                    A{i} = MPO_local_update_2sites( A{i}, (nn-1)*Q + bb, transpose( reshape( U_osc_int{bb}, [Nb(bb)^2 * Nb(bb+1)^2, Nb(bb)^2 * Nb(bb+1)^2]  ) ) );                                                          
                end               
                end
                % compress (faster than omitting and compress in system ev.)
                [A{i}, bonderr_tmp, ~] = MPO_compression_adap( A{i}, comp_tol , Nbe.^2 , chimax);        

                %chivals{i}(it,:) = chi_vec;
                %comp_number = comp_number +1;

                bonderr = bonderr + bonderr_tmp;

            end
            
    else
       
        bonderr = 0;
        
        for i=1:numel(Aid)             
            for nn = 1:N
            for bb = 1:Q-1    
                % Ex: Q=3, M=6, isite= (1)-2, (2)-3, (4)-5, (5)-6
                A{i} = MPO_local_update_2sites( A{i}, (nn-1)*Q + bb, transpose( reshape( U_osc_int{bb}, [Nb(bb)^2 * Nb(bb+1)^2, Nb(bb)^2 * Nb(bb+1)^2]  ) ) );                                                          
            end               
            end
            % compress (faster than omitting and compress in system ev.)
            [A{i}, bonderr_tmp, ~] = MPO_compression_adap( A{i}, comp_tol , Nbe.^2 , chimax);        

            %chivals{i}(it,:) = chi_vec;
            %comp_number = comp_number +1;

            bonderr = bonderr + bonderr_tmp;

        end
            
        
    end
        
end
