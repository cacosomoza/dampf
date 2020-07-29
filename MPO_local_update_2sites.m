% A_in:     input MPO 
% isite:    on which site the A-matrix changes (A_i^(isite))_jk
% y:        the local linear transformation given as a matrix. Size(Nb^4 X Nb^4) (enter as transpose)   

% There is an overall increase in the size of the tensors as a result of growing correlations between interacting oscillators.

function A_out = MPO_local_update_2sites(A_in, isite, y)
               
    D1  = size(A_in{isite}, 1);     % 1 or chi
    D2  = size(A_in{isite}, 2);     % 1 or chi
    d   = size(A_in{isite}, 3);     % Nb(isite)^2  
    D2p1 = size(A_in{isite+1}, 2);
    dp1  = size(A_in{isite+1}, 3);           
    A_out = A_in;
    
    % size of arrays before update
    %s1= numel(A_in{isite}); s2=numel(A_in{isite+1}); %fprintf( 'isite: %i ,numel(A_i{i}): %i, numel(A_i{i+1}): %i \n' , isite, s1, s2 )       
          
    X1 = reshape( permute( A_in{isite}, [1,3,2]), D1*d, [] );       
    X2 = reshape( A_in{isite+1}, D2,[] );  

    X = reshape( X1 * X2, [D1, d, D2p1, dp1] );        
    X = reshape( permute( X , [1,3,2,4] ), D1*D2p1, d*dp1 );
    
    X = X * y;
    
    Xtransf = reshape( permute( reshape(X, D1, D2p1, d, dp1 ),[1,3,2,4] ), D1*d, D2p1*dp1 );
        
    [Q,R] = qr( Xtransf, 0);       
               
    A_out{isite}   = permute( reshape( Q, D1, d,[] ), [1,3,2]);
    A_out{isite+1} = reshape( R, size(R,1), [],dp1 );
    
    % size of arrays after update (all larger except isite=1 and isite=5)
    %s1p=numel(A_out{isite}); s2p= numel(A_out{isite+1}); %fprintf( 'isite: %i ,numel(A_o{i}): %i, numel(A_o{i+1}): %i \n' ,isite, s1p, s2p )           
    %fprintf( 'site: %i, size change i: %i, size change i+1: %i \n' , isite, -s1+s1p, -s2+s2p )       
    
        
end

% % A_in:     input MPO 
% % isite:    on which site the A-matrix changes (A_i^(isite))_jk
% % y:        the local linear transformation given as a matrix. Size(Nb^4 X Nb^4)   
% 
% function A_out = MPO_local_update_2sites(A_in, isite, y)
%     
%     % A_in{isite}: 1 X chi X M^2, chi X chi X M^2 or chi X 1 X M^2
%     D1 = size(A_in{isite}, 1); % 1 or chi
%     D2 = size(A_in{isite}, 2); % 1 or chi
%     d  = size(A_in{isite}, 3); % Nb^2
%     
%     
% 
%     A_out = A_in;
%     
%     X = reshape( permute( reshape(reshape( permute( A_in{isite}, [1,3,2]),D1*d,[] ) * reshape( A_in{isite+1}, D2,[] ), D1,d,[],d ), [1,3,2,4] ), [],d^2 ) * y;
%     
%     [A1,A2] = qr( reshape( permute( reshape(X,D1,[],d,d ),[1,3,2,4] ), D1*d,[] ), 0);
%     
%     A_out{isite}   = permute( reshape( A1, D1, d,[] ), [1,3,2]);
%     A_out{isite+1} = reshape( A2, size(A2,1), [],d );
%     
% end