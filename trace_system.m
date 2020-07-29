    
function rhoscs = trace_system( A, Aid, it, tracefactor, rhoscs,  N, Q, oide, Nb, Nbe, getreducedosc, pop_idx, basis, parallel )


        if Nbe > 1
                % Total Reduced vibrational state (trace out electronic dof)
                rhov=A{1};
                for i=2:N 
                    rhov = MPO_sum( rhov, A{Aid==pop_idx(i)}, Nbe.^2 );
                end

                % Single oscillator reduced states
                for oo=1:numel(getreducedosc)

                    oscid = getreducedosc(oo);          % total index for selected oscillators to monitor            
                    coeffs = zeros(Nbe(oscid)^2,1);     % coefficients of selected rho_v_o in its own basis

                    % left osc            
                    if oscid == 1 

                        rest=rhov{2}(:,:,1);
                        for i=3:N*Q 
                            rest = rest * rhov{i}(:,:,1);
                        end
                        for ll=1:Nbe(oscid)^2  
                            coeffs(ll)  = rhov{1}(:,:,ll) * rest;
                        end

                    % right osc    
                    elseif oscid == N*Q 

                        rest = rhov{1}(:,:,1);
                        for i=2:N*Q-1 
                            rest = rest * rhov{i}(:,:,1);
                        end
                        for ll=1:Nbe(oscid)^2
                           coeffs(ll) = rest * rhov{N*Q}(:,:,ll); 
                        end

                    % middles oscs.
                    else

                        restL = rhov{1}(:,:,1);                
                        for i = 2:oscid-1                   
                           restL = restL * rhov{i}(:,:,1); 
                        end

                        restR = rhov{oscid+1}(:,:,1);
                        for i = oscid+1+1:N*Q 
                           restR = restR * rhov{i}(:,:,1); 
                        end                

                        for ll=1:Nbe(oscid)^2
                           coeffs(ll) = restL * rhov{oscid}(:,:,ll) * restR;
                        end

                    end                        
                    vibtracefactor = 1; % it is the product of trace(o{1}) of every oscillator except the one that remains
                    for qq=1:N*Q
                        if qq == oscid
                            % do nothing
                        else
                            vibtracefactor = vibtracefactor * trace( basis{oide(qq)}{1} ); %t The trace factor needs to consider that different modes have different basis and therefore diff o{1}
                        end
                    end                                                           
                    for ll=1:Nbe(oscid)^2            
                         rhoscs{oo}(it,:) = rhoscs{oo}(it,:) + reshape( vibtracefactor * coeffs(ll) * basis{oide(oscid)}{ll}, 1,[] );                  
                         % normalize state individual vibrational states
                          %whichosc=oide(getreducedosc(oo));
                          %%n_idx= diag( reshape(find(ones(Nb(whichosc),Nb(whichosc))),[Nb(whichosc),Nb(whichosc)]) );
                          %roscs{oo}(it,:) = rhoscs{oo}(it,:) / sum(rhoscs{oo}(it,n_idx),2);                 
                    end

                end                                   

                % Normalize total vibrational state for correct calculation of Mandel parameter
                vibtrace = rhov{1}(:,:,1);
                for i=2:N*Q
                    vibtrace = vibtrace * rhov{i}(:,:,1);
                end
                vibtrace = vibtrace * tracefactor;
                rhov = MPO_scalar(1/vibtrace, rhov);

        end                               

end