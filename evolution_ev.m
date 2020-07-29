
function A = evolution_ev( A, together, Aid, N, Q, Nb, Nbe, DDswch, DDleft, DDrght, parallel )

    if parallel
    
                if together % w1 w2 w1 w2                 

                        A_temp = cell(N,N);
                        parfor i = 1:numel(Aid)

                            [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables

                            for qq=1:Q
                                % populations
                                if j == k
                                    if qq==1                    
                                        A_temp{i} = MPO_local_update( A{i}     , Q*(j-1) + qq, DDswch{qq} );            
                                    else                    
                                        A_temp{i} = MPO_local_update( A_temp{i}, Q*(j-1) + qq, DDswch{qq} );            
                                    end
                                % coherences    
                                else
                                    if qq==1                    
                                        A_temp{i} = MPO_local_update( MPO_local_update( A{i}    ,  Q*(j-1) + qq, DDleft{qq} ), Q*(k-1) + qq, DDrght{qq} );            
                                    else                    
                                        A_temp{i} = MPO_local_update( MPO_local_update( A_temp{i}, Q*(j-1) + qq, DDleft{qq} ), Q*(k-1) + qq, DDrght{qq} );            
                                    end
                                end
                            end
                        end        
                        A = A_temp;    

                    else % w1 w1 w2 w2

                        A_temp = cell(numel(Aid),1);
                        parfor i = 1:numel(Aid)            

                            [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables

                            for qq=1:Q
                                if j == k
                                    if qq==1
                                        % populations
                                        A_temp{i} = MPO_local_update(A{i}, (qq-1)*N + j, DDswch{qq});            
                                    else
                                        % populations
                                        A_temp{i} = MPO_local_update(A_temp{i}, (qq-1)*N + j, DDswch{qq});            
                                    end
                                else
                                    if qq==1
                                        % coherences
                                        A_temp{i} = MPO_local_update( MPO_local_update( A{i}, (qq-1)*N + j, DDleft{qq} ), (qq-1)*N + k, DDrght{qq} );            
                                    else
                                        % coherences
                                        A_temp{i} = MPO_local_update( MPO_local_update( A_temp{i}, (qq-1)*N + j, DDleft{qq} ), (qq-1)*N + k, DDrght{qq} );            
                                    end
                                end
                            end
                        end      
                        A = A_temp; 
                end
                
    else

                 if together % w1 w2 w1 w2                 

                        A_temp = cell(N,N);
                        for i = 1:numel(Aid)

                            [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables

                            for qq=1:Q
                                % populations
                                if j == k
                                    if qq==1                    
                                        A_temp{i} = MPO_local_update( A{i}     , Q*(j-1) + qq, DDswch{qq} );            
                                    else                    
                                        A_temp{i} = MPO_local_update( A_temp{i}, Q*(j-1) + qq, DDswch{qq} );            
                                    end
                                % coherences    
                                else
                                    if qq==1                    
                                        A_temp{i} = MPO_local_update( MPO_local_update( A{i}    ,  Q*(j-1) + qq, DDleft{qq} ), Q*(k-1) + qq, DDrght{qq} );            
                                    else                    
                                        A_temp{i} = MPO_local_update( MPO_local_update( A_temp{i}, Q*(j-1) + qq, DDleft{qq} ), Q*(k-1) + qq, DDrght{qq} );            
                                    end
                                end
                            end
                        end        
                        A = A_temp;    

                    else % w1 w1 w2 w2

                        A_temp = cell(numel(Aid),1);
                        for i = 1:numel(Aid)            

                            [j,k] = ind2sub([N,N],Aid(i)); % temporary vaariables

                            for qq=1:Q
                                if j == k
                                    if qq==1
                                        % populations
                                        A_temp{i} = MPO_local_update(A{i}, (qq-1)*N + j, DDswch{qq});            
                                    else
                                        % populations
                                        A_temp{i} = MPO_local_update(A_temp{i}, (qq-1)*N + j, DDswch{qq});            
                                    end
                                else
                                    if qq==1
                                        % coherences
                                        A_temp{i} = MPO_local_update( MPO_local_update( A{i}, (qq-1)*N + j, DDleft{qq} ), (qq-1)*N + k, DDrght{qq} );            
                                    else
                                        % coherences
                                        A_temp{i} = MPO_local_update( MPO_local_update( A_temp{i}, (qq-1)*N + j, DDleft{qq} ), (qq-1)*N + k, DDrght{qq} );            
                                    end
                                end
                            end
                        end      
                        A = A_temp; 
                 end

    end
        
end
