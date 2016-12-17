function [E, serValSum] = stateEvolution(B, R, noiseParam, n, channel, monteCarloIterationCount)
    % state evolution for non-spatially coupled SS codes

    % noiseParam is either SNR or flip probability 
    if strcmpi(channel, 'bsc')
        Q_p = @(p,E) exp(-p.^2./(2.*E))./sqrt(2*pi.*E); % derivative of Q-function wrt p
        fisher = @(p,E) (Q_p(p,E) - 2*noiseParam.*Q_p(p,E)).^2 ./...
        ((qfunc(-p./sqrt(E)) + noiseParam - 2*noiseParam.*qfunc(-p./sqrt(E))).*(1-(qfunc(-p./sqrt(E)) + noiseParam - 2*noiseParam.*qfunc(-p./sqrt(E)))));
        fisher_int = @(p,E) fisher(p,E) .* exp(-p.^2./(2.*(1-E)))./sqrt(2*pi.*(1-E));
        fisher_exp = @(E) integral(@(p)fisher_int(p,E),-Inf,Inf,'ArrayValued',true); % expectation of fisher info
    elseif strcmpi(channel, 'bec')
        Q = @(p,E) qfunc(-p/sqrt(E));
        Q_p = @(p,E) exp(-(p^2)/(2*E))/sqrt(2*pi*E); % derivative of Q-function wrt p
        fisher = @(p,E) ((Q_p(p,E))^2)*(1-noiseParam)/max(Q(p,E)-Q(p,E)^2,1e-20);
        fisher_int = @(p,E) fisher(p,E) * exp(-(p^2)/(2*(1-E)))/sqrt(2*pi*(1-E));
        fisher_exp = @(E) integral(@(p) fisher_int(p,E),30*E-30,-30*E+30,'ArrayValued',true); % expectation of fisher info
    end


    E = zeros(n+1, 1);
    serValSum = zeros(n+1, 1);
    E(1) = 1;
    serValSum(1)=1;
    for t=1:n
        switch channel
            case 'awgnc'
                sigma = sqrt(R*(E(t) + 1/noiseParam));
            case 'bec'
                sigma = sqrt(R/fisher_exp(E(t)));
                if E(t) < 6e-5
                    E(t+1:end) = 0;
                    break;
                elseif isnan(sigma)
                    sigma = sqrt(R/(5*(1-noiseParam)/(pi*sqrt(2*pi))));
                end
            case 'bsc'
                sigma = sqrt(R/fisher_exp(E(t)));
                if E(t) < 6e-5
                    E(t+1:end) = 0;
                    break;
                elseif E(t) == 1
                    sigma = sqrt(R/fisher_exp(0.9999999));
                end
        end
        
        E(t+1) = 0;
        serValSum(t+1) = 0;
        for i = 1:monteCarloIterationCount
            z = randn(B,1);
            f_1 = 1/(1+exp(-log2(B)/(sigma^2))*sum(exp((z(2:end)-z(1))*sqrt(log2(B))/sigma))); % the first bit, which is supposed to be decoded as 1
            f_2 = zeros(B-1,1);
            for j=2:B
                f_2(j-1) = 1/(1+exp(log2(B)/(sigma^2)+(z(1)-z(j))*sqrt(log2(B))/sigma)+sum(exp((z([2:j-1 j+1:end])-z(j))*sqrt(log2(B))/sigma))); % the other bits, which are supposed to be decoded as all 0s
            end
            E(t+1) = E(t+1) + ((f_1-1)^2 +sum(f_2.^2)); % mse calculation
            serValSum(t+1) = serValSum(t+1) + double(max(f_2) > f_1); % ser calculation
        end
        E(t+1) = E(t+1) / monteCarloIterationCount;
        serValSum(t+1) = serValSum(t+1) / monteCarloIterationCount;
    end
end
