function Fu = potential(B, R, channel, noiseParam, E, monteCarloIterationCount)

    switch channel
        case 'awgnc'
            sigma = sqrt(R*(E + 1/noiseParam));
            secondExpectation = (-log(E*noiseParam+1)-log(noiseParam) + log(pi) + log(2) + 1)/(2*log(2));
            Uu = -E/(2*log(2)*sigma^2)-secondExpectation/R;
        case 'bsc'
            Q = @(p,E) qfunc(-p/sqrt(E));
            Q_p = @(p,E) exp(-(p^2)/(2*E))/sqrt(2*pi*E); % derivative of Q-function wrt p
            fisher = @(p,E) ((Q_p(p,E) - 2*noiseParam*Q_p(p,E))^2) / ((Q(p,E) + noiseParam - 2*noiseParam*Q(p,E))*(1 - Q(p,E) - noiseParam + 2*noiseParam*Q(p,E)));
            fisher_int = @(p,E) fisher(p,E) * exp(-(p^2)/(2*(1-E)))/sqrt(2*pi*(1-E));
            fisher_exp = @(E) integral(@(p) fisher_int(p,E),30*E-30,-30*E+30,'ArrayValued',true); % expectation of fisher info
            sigma = sqrt(R/fisher_exp(E));
            secondExpectation = integral(@(z) ...
                ((noiseParam*qfunc(-z*sqrt(1-E)/sqrt(E)) + (1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E))) * ...
                log2(noiseParam*qfunc(-z*sqrt(1-E)/sqrt(E)) + (1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E))) + ...
                ((1-noiseParam)*qfunc(-z*sqrt(1-E)/sqrt(E)) + noiseParam*qfunc(z*sqrt(1-E)/sqrt(E))) * ...
                log2((1-noiseParam)*qfunc(-z*sqrt(1-E)/sqrt(E)) + noiseParam*qfunc(z*sqrt(1-E)/sqrt(E)))) * ...
                exp(-(z^2)/2)/sqrt(2*pi), ...
                -30, 30,'ArrayValued', true);
            Uu = -E/(2*log(2)*sigma^2)-secondExpectation/R;
        case 'bec'
            Q = @(p,E) qfunc(-p/sqrt(E));
            Q_p = @(p,E) exp(-(p^2)/(2*E))/sqrt(2*pi*E); % derivative of Q-function wrt p
            fisher = @(p,E) ((Q_p(p,E))^2)*(1-noiseParam)/(max(Q(p,E),1e-18)*max((1-Q(p,E)),1e-18));
            fisher_int = @(p,E) fisher(p,E) * exp(-(p^2)/(2*(1-E)))/sqrt(2*pi*(1-E));
            fisher_exp = @(E) integral(@(p) fisher_int(p,E),30*E-30,-30*E+30,'ArrayValued',true); % expectation of fisher info
            sigma = sqrt(R/fisher_exp(E));
            secondExpectation = integral(@(z) ...
                (noiseParam*log2(noiseParam) + ...
                (1-noiseParam)*qfunc(-z*sqrt(1-E)/sqrt(E))*log2(max(1e-300, (1-noiseParam)*qfunc(-z*sqrt(1-E)/sqrt(E)))) + ...
                (1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E))*log2(max(1e-300, (1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E))))) * ...
                exp(-(z^2)/2)/sqrt(2*pi), ...
                -30, 30,'ArrayValued', true);
            Uu = -E/(2*log(2)*sigma^2)-secondExpectation/R;
		case 'zc'
			Q = @(p,E) qfunc(-p/sqrt(E));
			Q_p = @(p,E) exp(-(p^2)/(2*E))/sqrt(2*pi*E); % derivative of Q-function wrt p
			fisher = @(p,E) ((Q_p(p,E)*(1-noiseParam))^2)/max(Q(p,E)+noiseParam*(1-Q(p,E)),1e-18) + ...
				(Q_p(p,E)^2)*(1-noiseParam)/max(1-Q(p,E), 1e-18);
			fisher_int = @(p,E) fisher(p,E) * exp(-(p^2)/(2*(1-E)))/sqrt(2*pi*(1-E));
			fisher_exp = @(E) integral(@(p) fisher_int(p,E),30*E-30,-30*E+30,'ArrayValued',true); % expectation of fisher info
			sigma = sqrt(R/fisher_exp(E));
            secondExpectation = integral(@(z) ...
                ((1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E))*log2(max(1e-300, (1-noiseParam)*qfunc(z*sqrt(1-E)/sqrt(E)))) + ...
                (noiseParam*qfunc(z*sqrt(1-E)/sqrt(E))+qfunc(-z*sqrt(1-E)/sqrt(E))) * ...
                log2(max(1e-300, noiseParam*qfunc(z*sqrt(1-E)/sqrt(E))+qfunc(-z*sqrt(1-E)/sqrt(E))))) * ...
                exp(-(z^2)/2)/sqrt(2*pi), ...
                -30, 30,'ArrayValued', true);
            Uu = -E/(2*log(2)*sigma^2)-secondExpectation/R;
    end
            
    switch B
        case 2
            Su = log2(exp(1))*(- 1 - log(2) + integral(@(z) log(1+exp(-1/(sigma^2)-z*sqrt(2)/sigma))*exp(-(z^2)/2)/sqrt(2*pi),-30,30,'ArrayValued',true));
        otherwise
            s = zeros(1,B);
            s(1) = 1;
            monteCarloSum = 0;
            for i = 1:monteCarloIterationCount
                z = randn(1,B);
                insideLog = 0;
                for k=1:B
                    x = zeros(1,B);
                    x(k) = 1;
                    Q = exp(-log2(B)*(norm(x-(s+z*sigma/sqrt(log2(B))))^2)/(2*sigma^2));
                    insideLog = insideLog + Q/B;
                end
                monteCarloSum = monteCarloSum + log(insideLog)/log(B);
            end
            Su = monteCarloSum / monteCarloIterationCount;
    end
    Fu = Uu - Su;
end
