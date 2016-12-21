% GAMP algorithm for Hadamard matrices, stores the results of all iterations

% noiseParam is either SNR or erasure/flip probability
shat = 0;
xhat = zeros(N, n+1);
Tx = ones(N, 1).*(1/B);

%% GAMP
for t=1:n
    %% output linear step
    Tp = MultSeededHadamardSquarred(Tx, J, Lr, Lc, subsectionRowSizes, N/Lc);
    phat = MultSeededHadamard(xhat(:,t), J, Lr, Lc, subsectionRowSizes, N/Lc, rp, []) - Tp.*shat;
    
    %% output nonlinear step
    switch channel
        case 'awgnc'
            shat = (y-phat)./(Tp+1/noiseParam);
            Ts = 1./(Tp+1/noiseParam);
        case 'bsc'
            Z = erfc(phat./sqrt(2*Tp)).*((1-noiseParam).*double(y==-1)+noiseParam.*double(y==1))/2 + ...
                (1+erf(phat./sqrt(2*Tp))).*((1-noiseParam).*double(y==1)+noiseParam.*double(y==-1))/2;
            zhat = (-exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)./sqrt(2*pi) + phat.*erfc(phat./sqrt(2*Tp))/2).*...
                ((1-noiseParam).*double(y==-1)+noiseParam.*double(y==1))+...
                (exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)./sqrt(2*pi) + phat.*(1/2+erf(phat./sqrt(2*Tp))/2)).*...
                ((1-noiseParam).*double(y==1)+noiseParam.*double(y==-1));
            zhat = zhat./Z;
            shat = (zhat - phat)./Tp;
            Ezsquare = (-exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)./sqrt(2*pi) + (phat.^2+Tp).*erfc(phat./sqrt(2*Tp))/2).*...
                ((1-noiseParam)*double(y==-1)+noiseParam*double(y==1)) +...
                (exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)./sqrt(2*pi) + (phat.^2+Tp).*(1/2+erf(phat./sqrt(2*Tp))/2)).*...
                ((1-noiseParam)*double(y==1)+noiseParam*double(y==-1));
            Ezsquare = Ezsquare./Z;
            Ts = (Tp - Ezsquare + zhat.^2)./(Tp.^2);
        case 'bec'
            Z = erfc(phat./sqrt(2*Tp)).*((1-noiseParam).*double(y==-1))/2 + ...
                (1+erf(phat./sqrt(2*Tp))).*((1-noiseParam).*double(y==1))/2 + ...
                noiseParam.*double(y==0);
            zhat = (-exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)./sqrt(2*pi) + phat.*erfc(phat./sqrt(2*Tp))/2).*...
                ((1-noiseParam).*double(y==-1))+...
                (exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)./sqrt(2*pi) + phat.*(1/2+erf(phat./sqrt(2*Tp))/2)).*...
                ((1-noiseParam).*double(y==1))+...
                noiseParam.*double(y==0).*phat;
            zhat = zhat./Z;
            shat = (zhat - phat)./Tp;
            Ezsquare = (-exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)./sqrt(2*pi) + (phat.^2+Tp).*erfc(phat./sqrt(2*Tp))/2).*...
                ((1-noiseParam)*double(y==-1)) +...
                (exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)./sqrt(2*pi) + (phat.^2+Tp).*(1/2+erf(phat./sqrt(2*Tp))/2)).*...
                ((1-noiseParam)*double(y==1)) +...
                noiseParam.*double(y==0).*(Tp+phat.^2);
            Ezsquare = Ezsquare./Z;
            Ts = (Tp - Ezsquare + zhat.^2)./(Tp.^2);
        case 'zc'
            Z = erfc(phat./sqrt(2*Tp)).*((1-noiseParam).*double(y==-1) + noiseParam.*double(y==1))/2 + ...
                (1+erf(phat./sqrt(2*Tp))).*double(y==1)/2;
            zhat = (-exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)/sqrt(2*pi) + phat.*erfc(phat./sqrt(2*Tp))/2) .* ...
                    ((1-noiseParam).*double(y==-1) + noiseParam.*double(y==1)) +  ...
                   (exp(-(phat.^2)./(2*Tp)).*sqrt(Tp)/sqrt(2*pi) + phat.*(1+erf(phat./sqrt(2*Tp)))/2) .* ...
                    double(y==1);
            zhat = zhat ./ Z;
            shat = (zhat - phat)./Tp;
            Ezsquare = (-exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)/sqrt(2*pi) + (phat.^2+Tp).*erfc(phat./sqrt(2*Tp))/2) .* ...
                        ((1-noiseParam).*double(y==-1) + noiseParam.*double(y==1)) + ...
                       (exp(-(phat.^2)./(2*Tp)).*phat.*sqrt(Tp)/sqrt(2*pi) + (phat.^2+Tp).*(1+erf(phat./sqrt(2*Tp)))/2) .* ...
                        double(y==1);
            Ezsquare = Ezsquare./Z;
            Ts = (Tp - Ezsquare + zhat.^2)./(Tp.^2);
    end
    
    %% input linear step
    Tr = max(1./(MultSeededHadamardTransposeSquarred(Ts, J, Lr, Lc, subsectionRowSizes, N/Lc)), 1e-100); % min value is set due to the numerical problems of BEC and BSC
    rhat = xhat(:,t) + Tr.*(MultSeededHadamardTranspose(shat, J, Lr, Lc, subsectionRowSizes, N/Lc, rp, []));
    
    %% input nonlinear step
    xhat(:,t+1) = gin(rhat, Tr, B);
    % if the resulting sections contain NaN values, we just reload the previous results.
    for secNo=1:L
        if isnan(sum(xhat((secNo-1)*B+1:secNo*B, t+1)))
            xhat((secNo-1)*B+1:secNo*B,t+1) = xhat((secNo-1)*B+1:secNo*B,t);
        end
    end
    Tx = xhat(:,t+1)-xhat(:,t+1).^2;
end