function x = hardDecision(xhat, B)
    % maps xhat to SS codex such that the largest value in each section is
    % set to 1, and the rest to 0. And does this for each iteration result.

    x = zeros(size(xhat));
    for iter = 1:size(xhat,2)
        for j=1:B:size(xhat,1)
            x(j:j+B-1,iter) = double(xhat(j:j+B-1,iter)==max(xhat(j:j+B-1,iter)));
        end
    end
end