function y = ser(c1, c2, B)
    L = length(c1)/B;
    correct = 0;
    for i=1:B:length(c1)
        correct = correct + double(max(abs(c1(i:i+B-1)-c2(i:i+B-1)))<0.1);
    end
    y = 1 - correct / L;
end