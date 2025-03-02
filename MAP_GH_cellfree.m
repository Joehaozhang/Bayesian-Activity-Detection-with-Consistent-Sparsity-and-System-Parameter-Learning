function [G_hat, z] = MAP_GH_cellfree(Y, S)
% GHVI
% Variational Inference Activity Detection and Channel Estimation
% Author: Hao Zhang
%
% Last updated: 2025/03/01
% Usage:
% [G_hat, z] = VIAD(Y, S)
% Solves for X = S{Gamma}H in
% Y = X + N = S{Gamma}H + N
% N is Gaussian noise
% -----------------------------------------------------------------------
%                                 INPUTS                                |
% Y: input matrix (received signals in BS)                              |
% S: pilot signals transmitted                                          |
% verbose: output the progress? (0/1) default: 1.                       |
% MAXITER: max number of iterations. default: 200.                      |
% Threshold: algorithm stop threshold. default: 1e-4.                   |
% -----------------------------------------------------------------------
%                                 OUTPUTS                               |
% X: the reconsturcted signals                                          |
% G: channel matrix where each row represents the corresponding channel |
% vector                                                                |
% -----------------------------------------------------------------------
%% Initialization 
[L, M, K] = size(Y); % L: sequence length, M: antenna, K: BS
LM = L*M;
[~, N] = size(S);
Y2sum  = trace(Y(:,:,1) * Y(:,:,1)');
S2sum  = trace(S * S');
scale2 = Y2sum / (S2sum * LM);

% Random initialization for Gamma and H
Gamma = rand(N,K) * 1; 
H = repmat(zeros(N,M),[1 1 K]);
for k=1:K
    H(:,:,k) = complex(randn(N,M),randn(N,M));
end
G_hat = repmat(zeros(N,M),[1 1 K]);

% Calculate X
X = repmat(zeros(L,M),[1 1 K]);
for k=1:K
    X(:,:,k) = S * diag(Gamma(:,k)) * H(:,:,k);
end
%% Variational distribution initialization
% GH distribution parameters
z_inv   = K*ones(N,1);
z       = 1 ./ z_inv;
a        = 8e-4*ones(N,1);
% b        = 10*rand(N,1);
b       = 10*ones(N,1);

% Gamma noise distribution
tau    = 1./scale2;

% Gamma noise prior
e_0      = 1e-6;
f_0      = 1e-6;

% Algorithm Parameter
error = zeros(200); % Error record
MAXITER     = 200;  % Maximal iterations
verbose     = 1;   % Progress display
Threshold = 1e-4;% Algorithm break threshold
%% Iterations
tic
for it = 1:MAXITER
    %% Save X
    old_X = X;
    %% Update Gamma
    for k=1:K
        YminusX = Y(:,:,k) - old_X(:,:,k);
        Y_minusi = YminusX;
        for i=1:N
            HiT = H(i,:,k);
            Si = S(:,i);
            Y_minusi =  Y_minusi + Gamma(i,k) * Si * HiT;
%             Sigma_L(:,:,k,i) = ((beta * (Si' * Si) * (HiT * HiT') + z_inv(i))).^(-1);
            Gamma(i,k) = tau * ((tau * (Si' * Si) * (HiT * HiT') + z_inv(i)/2)).^(-1) ...
                * trace(real(HiT' * Si' * Y_minusi)) ;
            Y_minusi = Y_minusi - Gamma(i,k) * Si * HiT;
        end
    end
    %% Update H
    for k=1:K
        X(:,:,k) = S * diag(Gamma(:,k)) * H(:,:,k);
    end
    for k=1:K
        YminusXT = (Y(:,:,k) - X(:,:,k)).';
%         YT_minusi = YminusXT;
        for i=1:N
            Si = S(:,i);
            Gamma_i = Gamma(i,k);
            YT_minusi = YminusXT + (Gamma_i * Si * H(i,:,k)).';
            H(i,:,k) = tau * Gamma_i * (tau * (Gamma_i.^2) * (Si' * Si) + 1).^(-1) * eye(M) * YT_minusi * conj(Si);
        end
    end
    %% Update X
    for k=1:K
        X(:,:,k) = S * diag(Gamma(:,k)) * H(:,:,k);
    end
    %% Update z
    for i=1:N
        z(i) = ((lambda_0-1-K/2)+sqrt((lambda_0-1-K/2)^2+a(i)*b(i)+2*a(i)*sum(Gamma(i,:).^2)))/a(i);
        z_inv(i) = 1./z(i);
    end
    z     = real(z);
    z_inv = real(z_inv);
    %% Update a
    for n=1:N
        a_old = a(n);
        for i=1:100
            DfDa = (1e-6 - 2)/(2*a(n)) - 1e-6 - 1/2*z(n) + (5*besselk(-5, (a(n)*b(n))^(1/2))/(2*a(n)) + (b(n)*besselk(-4, (a(n)*b(n))^(1/2)))/(2*(a(n)*b(n))^(1/2)))/besselk(-5, (a(n)*b(n))^(1/2));
            a(n) = real(a(n) + 0.05 * DfDa);
            if abs(a(n)-a_old)<0.00001
                break;
            end
            a_old = a(n);
        end
    end
    %% Update beta
    err = 0;
    for k=1:K 
        delta = Y(:,:,k) - X(:,:,k);
        err = err + trace(delta' * delta);
    end
    err_novariance = err;
    error(it) = err;
    tau = (K*LM + e_0)./(err_novariance + f_0);
    %% Display progress
    relative_change = 0;
    for k=1:K
        delta_X = old_X(:,:,k) - X(:,:,k);
        relative_change = relative_change + trace(delta_X' * delta_X) / (LM);
    end
    if verbose
        % For synthetic simulations
        recons_err = 0;
        for k=1:K
            delta_XY = X(:,:,k) - Y(:,:,k);
            recons_err = recons_err + trace(delta_XY' * delta_XY)...
                /trace(Y(:,:,k)' * Y(:,:,k));
        end
        recons_err = recons_err/K;
        fprintf('Method: MAP, iter %d: error = %g, tau = %g, relative change = %g\n', it, recons_err, tau, relative_change);
    end
    %% Check convergence
    if it > 20 && (relative_change < Threshold || error(it)>error(it-1))%
        break;
    end
end
for k=1:K
    G_hat(:,:,k) = diag(Gamma(:,k)) * H(:,:,k);
end
end