function expectation = gin(r, Tr, B)

    % g_in function for GAMP's input nonlinear step

    L = length(r)/B;
    numerator = exp((2*r-1)./(2*Tr));
    denominator = zeros(B*L,1);
    for i=1:B:B*L
        denominator(i:i+B-1) = sum(numerator(i:i+B-1));
    end
    expectation = numerator ./ denominator;

end