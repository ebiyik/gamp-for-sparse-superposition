function code = sparseSuperposition(B,L)
    % creates a random SS code with the specified B and L values
    for i=1:L
        temp = rand(B,1);
        code((i-1)*B+1:i*B) = (temp==max(temp));
    end
    code = double(code');
end