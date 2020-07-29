function A = MPO_conjugate(A, s, Nbe2, oide) 

% Input MPO A (cell of N 3d matrices) 
% s, signs of M^2 operators, according to their symmetry

for j = 1:length(A)
    
    A{j} = conj(A{j});
    
    for k = 1:Nbe2(j)
        if k > 1
        A{j}(:,:,k) = (2*s{oide(j)}(k)-1)*A{j}(:,:,k);
        else
            %A{j}(:,:) = (2*s{oide(j)}(k)-1)*A{j}(:,:);
        end
    end
end
