function A = freebath_damping( A ,together, Aid,  N, Q, damp, oide, BBath, Damp_down, Damp_up, parallel )
    
    if parallel

        parfor i = 1:numel(Aid)               
            for l = 1:N*Q % Mode 
                % Free Bath Evolution
                A{i} = MPO_local_update( A{i}, l, BBath{oide(l)} );                  
                % Oscillator Damping
                if damp(oide(l))==1                
                    if together % surrogate oscs                
                        A{i} = MPO_local_update( A{i}, l, Damp_down{oide(l)} );                  
                    else
                        A{i} = MPO_local_update( A{i}, l, Damp_up{oide(l)}   );
                        A{i} = MPO_local_update( A{i}, l, Damp_down{oide(l)} );                                  
                    end
                end
            end              
        end
     
    else
        
            for i = 1:numel(Aid)

                for l = 1:N*Q % Mode 
                    % Free Bath Evolution
                    A{i} = MPO_local_update( A{i}, l, BBath{oide(l)} );                  
                    % Oscillator Damping
                    if damp(oide(l))==1                
                        if together % surrogate oscs                
                            A{i} = MPO_local_update( A{i}, l, Damp_down{oide(l)} );                  
                        else
                            A{i} = MPO_local_update( A{i}, l, Damp_up{oide(l)}   );
                            A{i} = MPO_local_update( A{i}, l, Damp_down{oide(l)} );                                  
                        end
                    end                
                end                                  
            end                 
    end
          
end