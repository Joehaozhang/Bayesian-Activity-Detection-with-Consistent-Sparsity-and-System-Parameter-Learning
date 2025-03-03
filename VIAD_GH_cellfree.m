function [G_hat,z] = VIAD_GH_cellfree(Y, S)
% GHVI
% Variational Inference Activity Detection and Channel Estimation
% Author: Hao Zhang
%
% Last updated: 2022/07/08
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
scale  = sqrt(scale2);

% Random initialization for Gamma and H
Gamma = randn(N,K) * scale;
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

% Initialize variance
Sigma_L = repmat( scale * eye(1), [1 1 K N]);
Sigma_H = repmat( scale * eye(M,M), [1 1 K N]);
%% Variational distribution initialization
% GH distribution parameters
z_inv   = K * ones(N,1);
z       = 1 ./ z_inv;
a        = zeros(N,1);
b        = zeros(N,1);
lambda   = zeros(N,1);

% GIG distribution prior
a_0      = abs(ones(N,1));
b_0      = 1e-6;
lambda_0 = 1e-6;
kappa_a1 = -lambda_0/2 + 1e-6;
kappa_a2 = 1e-6;

% Gamma noise distribution
tau    = 1./scale2;

% Gamma noise prior
e_0      = 1e-6;
f_0      = 1e-6;

% Algorithm Parameter
error = zeros(200); % Error record
MAXITER     = 100;  % Maximal iterations
verbose     = 1;   % Progress display
Threshold = 1e-4;% Algorithm break threshold
%% Iterations
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
            Y_minusi = Y_minusi + Gamma(i,k) * Si * HiT;
            Sigma_L(:,:,k,i) = (2*(tau * (Si' * Si) * (HiT * HiT' ...
                + trace(Sigma_H(:,:,k,i))) + z_inv(i)/2)).^(-1);
            Gamma(i,k) = 2*tau * Sigma_L(:,:,k,i) ...
                * trace(real(HiT' * Si' * Y_minusi));
            Y_minusi = Y_minusi - Gamma(i,k) * Si * HiT;
        end
    end
    %% Update H
    for k=1:K
        X(:,:,k) = S * diag(Gamma(:,k)) * H(:,:,k);
    end
    for k=1:K
        YminusXT = (Y(:,:,k) - X(:,:,k)).';
        YT_minusi = YminusXT;
        for i=1:N
            Si = S(:,i);
            Gamma_i = Gamma(i,k);
            YT_minusi = YT_minusi + (Gamma_i * Si * H(i,:,k)).';
            Sigma_H(:,:,k,i) = (tau * (Gamma_i.^2 + Sigma_L(:,:,k,i)) ...
                * (Si' * Si) ...
                + 1).^-1 * eye(M);
            H(i,:,k) = tau * Gamma_i * Sigma_H(:,:,k,i) * YT_minusi * conj(Si);
            YT_minusi = YT_minusi - (Gamma_i * Si * H(i,:,k)).';
        end
    end
    %% Update X
    for k=1:K
        X(:,:,k) = S * diag(Gamma(:,k)) * H(:,:,k);
    end
    %% Update z
    for i=1:N
        a(i)  = a_0(i);
        Sig_L = 0;
        for k=1:K
            Sig_L = Sig_L + Sigma_L(:,:,k,i);
        end
        b(i)      = b_0 + (sum(Gamma(i,:).^2) + Sig_L);
        lambda(i) = lambda_0-K/2;
        z_inv(i)  = (b(i)/a(i))^(-1/2) * (besselk(lambda(i)+1,sqrt(a(i) ...
            * b(i)))/besselk(lambda(i),sqrt(a(i) * b(i))));
        z(i)      = (b(i)/a(i))^(1/2)  * (besselk(lambda(i)+1,sqrt(a(i) ...
            * b(i)))/besselk(lambda(i),sqrt(a(i) * b(i))));
    end
    %% Update beta
    err = 0;
    for k=1:K
        delta = Y(:,:,k) - X(:,:,k);
        err = err + trace(delta' * delta);
    end

    for k=1:K
        for i=1:N
            Si = S(:,i);
            Hi = H(i,:,k);
            SigmaLi = Sigma_L(:,:,k,i);
            SigmaHi = Sigma_H(:,:,k,i);
            err = err + (Gamma(i,k)^2 * trace(SigmaHi) ...
                + SigmaLi * (Hi * Hi') ...
                + trace(SigmaLi * SigmaHi)) ...
                * (Si' * Si);
        end
    end
    error(it) = err;
    tau = (K*LM + e_0)./(err + f_0);
    %% update a_0
    for i=1:N
        a_0(i) = (kappa_a1 + lambda_0/2)/(kappa_a2 + z(i)/2);
    end
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
        fprintf('Method: VIGH, iter %d: error = %g, tau = %g, relative change = %g\n', it, recons_err, tau, relative_change);
    end
    %% Check convergence
    if it > 20 && (relative_change < Threshold)
        break;
    end
end
%% Generate channel estimation
for k=1:K
    G_hat(:,:,k) = diag(Gamma(:,k)) * H(:,:,k);
end
end