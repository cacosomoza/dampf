    function A_out = MPO_sum(A_in1, A_in2,Nbe)

    assert(length(A_in1) == length(A_in2));
    N = length(A_in1);
    %d = size(A_in1{1},3);
    
    
    A_out = cell(1,N);  
    A_out{1} = zeros(1, size(A_in1{1},2)+size(A_in2{1},2), Nbe(1));
    A_out{1}(1, 1:size(A_in1{1},2), :) = A_in1{1};
    A_out{1}(1, size(A_in1{1},2)+(1:size(A_in2{1},2)), :) = A_in2{1};
               
    for j = 2:N-1    
        D11 = size(A_in1{j}, 1);
        D12 = size(A_in1{j}, 2);
        D21 = size(A_in2{j}, 1);
        D22 = size(A_in2{j}, 2);
        A_out{j} = zeros(D11+D21, D12+D22, Nbe(j));
        A_out{j}(1:D11, 1:D12, :) = A_in1{j};
        A_out{j}(D11+1:D11+D21, D12+1:D12+D22, :) = A_in2{j};
    end
        
        A_out{N} = zeros(size(A_in1{N},1)+size(A_in2{N},1),1,Nbe(N));
        A_out{N}(1:size(A_in1{N},1), 1, :) = A_in1{N};
        A_out{N}(size(A_in1{N},1)+(1:size(A_in2{N},1)), 1, :) = A_in2{N};

end