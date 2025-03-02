[N,M,K,monte]       = size(G_hat_ghvi);
G_hat_dominant_ghvi = zeros(N,M,monte);
G_hat_dominant_map  = zeros(N,M,monte);
G_real_dominant    = zeros(N,M,monte);
DominantAP          = zeros(N,monte);

for j=1:monte
    Gnorm2   = G_real(:,:,:,j).*conj(G_real(:,:,:,j));
    norm2sum = sum(Gnorm2,2);
    [~,DominantAP(:,j)] = max(norm2sum,[],3);
end

for j=1:monte
    for n=1:N
        G_hat_dominant_ghvi(n,:,j) = G_hat_ghvi(n,:,DominantAP(n,j),j);
        G_hat_dominant_map(n,:,j)  = G_hat_map(n,:,DominantAP(n,j),j);
        G_real_dominant(n,:,j)    = G_real(n,:,DominantAP(n,j),j); 
    end
end
Gnorm2sum_real     = sum(G_real_dominant .* conj(G_real_dominant),2);
Gnorm2sum_hat_ghvi = sum(G_hat_dominant_ghvi .* conj(G_hat_dominant_ghvi),2);
Gnorm2sum_hat_map  = sum(G_hat_dominant_map .* conj(G_hat_dominant_map),2);