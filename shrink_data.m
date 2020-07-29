
% shrinks data output

function shrink_data(pathtofiles, cases)

    % input 1: path 
    % input 2: array with id

    for i=1:numel(cases)
        if isfile( sprintf(pathtofiles,cases(i)) )
            try
                load( sprintf(pathtofiles,cases(i) ) , 'rho_out' ,'tmaxfs','dtfs','sysID',...
                              'AccTruncerr','Truncerr','oide','getreducedosc','Nb','rhoscs','occup2','occup','Nt','chivals','N','vibrations_info');
                
                %load( sprintf(pathtofiles,cases(i) ) , 'Aid','repeats','rho_out' ,'tmaxfs','dtfs','sysID',...
                %              'AccTruncerr','Truncerr','oide','getreducedosc','Nb','rhoscs','occup2','occup','Nt','chivals','N','vibrations_info');
                
                
                save( sprintf(pathtofiles,cases(i) ) );
  
            catch
                fprintf('error loading data_%i.mat, skipped.\n', cases(i));
            end
            
        else
            fprintf('data_%i is not a file \n',cases(i)) 
            
        end
        
    end
    
    disp('finished!')
                  
     
end


