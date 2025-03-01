clear all
clc
%% Notation
% ----------------------------------------                      
% |D       |Area size                    |
% |M       |Number of AP                 |
% |K       |Number of potential users    |
% |epsilon |Active probability           |
% |L       |Pilot sequence length        |
% |N       |Number of antennas           |
% |S       |Pilot sequence              |
% |gamma   |Active indicator and pathloss|
% |H       |Small-scale fading           |
% |Y       |Received signals             |
% ----------------------------------------
%% Simulation setting
% Area Size(km)
D = 3; 

% Number of BS                                                         
M = 12;

% Total Number of Users in Cell Network 
K = 500;

% Number of Antennas at each AP 
N = 40;

% Binomial Activity percentage: Activity pattern
epsilon = 25/K;

% Pilot Sequence Length 
L = 25;

% Noise Parameter : bandwidth = 1 MHz
sigma_sqr_dBm = -109;
sigma_sqr = 10^((sigma_sqr_dBm-30)/10); % Real noise variance
sigma_sqrN = 1; % Normalized Sigma2

% Transmit Power Constraint in W
TxPow = 200e-3; % 23dBm

monte = 500;
%% Generate Signature Sequence (multiplied by transmit power)
% S = (1/sqrt(2)*complex(randn(L,K),randn(L,K)));
S = (1/sqrt(2*L)*complex(randn(L,K),randn(L,K)));
S = S./vecnorm(S,2,1);

%%    Locations of BS and Users
    % AP=unifrnd(-D/2,D/2,M,2); % Uniformly Allocate the AP locations
    % 20 AP
    % AP = [-2*D/5,3*D/8;-D/5,3*D/8;0,3*D/8;D/5,3*D/8;2*D/5,3*D/8;...
    %     -2*D/5,D/8;-D/5,D/8;0,D/8;D/5,D/8;2*D/5,D/8;...
    %     -2*D/5,-D/8;-D/5,-D/8;0,-D/8;D/5,-D/8;2*D/5,-D/8;...
    %     -2*D/5,-3*D/8;-D/5,-3*D/8;0,-3*D/8;D/5,-3*D/8;2*D/5,-3*D/8;];
    %     16 AP
    % AP = 1/4 * [-D,D;-D/6,D;D/6,D;D,D;...
    %     -D,D/6;-D/6,D/6;D/6,D/6;D,D/6;...
    %     -D,-D/6;-D/6,-D/6;D/6,-D/6;D,-D/6;...
    %     -D,-D;-D/6,-D;D/6,-D;D,-D];
%     12AP
    AP = [-D/3,3*D/8;0,3*D/8;D/3,3*D/8;-D/3,D/8;0,D/8;D/3,D/8;-D/3,-D/8;0,-D/8;D/3,-D/8;-D/3,-3*D/8;0,-3*D/8;D/3,-3*D/8];
%     10AP
    % AP = [-D/4,0;-D/4,-D/4;0,-D/4;D/4,-D/4;D/4,0;D/4,D/4;0,D/4;-D/4,D/4;0,0;0,0];
%     9 AP
%     AP = [-D/4,0;-D/4,-D/4;0,-D/4;D/4,-D/4;D/4,0;D/4,D/4;0,D/4;-D/4,D/4;0,0];
%     8 AP
    % AP = [-D/4,0;-D/4,-D/4;0,-D/4;D/4,-D/4;D/4,0;D/4,D/4;0,D/4;-D/4,D/4];
%     6 AP
%     AP = [-D/6,-D/6;D/6,-D/6;D/6,D/6;-D/6,D/6;-D/3,0;D/3,0];
%     5 AP
%     AP = [-D/6,-D/6;D/6,-D/6;D/6,D/6;-D/6,D/6;0,0];
%     4 AP
    % AP = [-D/4,-D/4;D/4,-D/4;D/4,D/4;-D/4,D/4];
%     3 AP
%     AP = [-D/6,-D/6;0,D/6;D/6,-D/6];
%     2 AP
    % AP = [-D/4,0;D/4,0];
%     1AP
     % AP = [0,0];
%% Number of Monte-Carlo simulations
runtime_vigh = 0;
runtime_map = 0;
%% Variables for evaluation initialization
Y_real = zeros(L,N,M);
Y      = Y_real;
gamma       = zeros(K,monte,M); 
G_real      = repmat(zeros(K,N),[1 1 M monte]);
% t           = zeros(1,monte);
Active_List = zeros(K,monte);
gamma_real  = zeros(K,monte);
gamma_hat   = zeros(K,monte);
Beta        = zeros(K,M);
Beta_true   = zeros(K,M);
Beta2       = zeros(K,M);
lambda_gmmvamp = repmat(zeros(K,N),[1 1 M monte]);
%% Estimation Initialization
G_hat_CAMP     = repmat(zeros(K,N),[1 1 M monte]);
P_a_CAMP       = zeros(K,monte);
conv_CAMP      = zeros(1,monte);
G_hat_CVAMP     = repmat(zeros(K,N),[1 1 M monte]);
P_a_CVAMP       = zeros(K,monte);
conv_CVAMP      = zeros(1,monte);
G_hat_AMP      = repmat(zeros(K,N),[1 1 M monte]);
G_hat_vigh  = repmat(zeros(K,N),[1 1 M monte]);
% G_hat_vigh_v2  = repmat(zeros(K,N),[1 1 M monte]);
% G_hat_vigh_v3  = repmat(zeros(K,N),[1 1 M monte]);
% G_hat_map   = repmat(zeros(K,N),[1 1 M monte]);
G_hat_gmmvamp = repmat(zeros(K,N),[1 1 M monte]);
G_hat_rowAMP  = repmat(zeros(K,N),[1 1 M monte]);
% z_hat_vigh = zeros(K,monte);
% z_inv_hat_vigh = zeros(K,monte);
%% Generate received signals
for i = 1:1:monte
%     Random active devices (Binoulli)
%     Ac_list = binornd(1,epsilon,[K,1]);
%     Active_List(:,i) = Ac_list;

%     Random active devices (fixed number)
    Ac_idx = randperm(K);
    Ac_list = zeros(K,1);
    Ac_list(Ac_idx(1:K*epsilon)) = 1;
%     Active_List(:,i) = Ac_list;

%     Index of active devices
    idx = find(Ac_list);
%%     Uniformly allocate the UT locations in a DxD area
    UT=unifrnd(-D/2,D/2,K,2);
%     UT=unifrnd(-1.5,1.5,K,2);
%%     Uniformly allocate the UT locations in a circular area
%     UT_angle = rand(K,1) * 2 * pi;
%     UT_r     = rand(K,1) * D/2    ;
%     UT(:,1)  = UT_r .* cos(UT_angle);
%     UT(:,2)  = UT_r .* sin(UT_angle);
    for m=1:M
        %     Pathloss (distance from M BS)
        dist = sqrt(sum((UT - AP(m,:)).^2,2));
        %     Fixed distance
%             dist = 0.5 + 0.5 * rand(K,1);
        Beta(:,m) = 0.2 * 10.^((-128.1 - 37.6*log10(dist))/10);
        Beta_true(:,m) = 0.2 * 10.^((-128.1 - 37.6*log10(dist)+2*rand(K,1))/10);
        %     Pathloss for normalized noise (variance=1)
        Beta(:,m) = Beta(:,m)/sigma_sqr;
    end

    maxBeta           = max(Beta,[],2);
    maxBeta2          = maxBeta;
    idx_snr           = find(maxBeta > 4);
    maxBeta2(idx_snr) = 4;
    idx_snr2 = find(maxBeta < 4);
    maxBeta2(idx_snr2) = 4;
    PowControl        = maxBeta2./maxBeta;
    Beta2             = PowControl.*Beta_true/sigma_sqr;
    idx2              = find(maxBeta2>0);
    idx               = intersect(idx,idx2);
    Active_List(idx,i)=1;
    for m=1:M
        %     Pathloss (distance from M BS)
%         dist = sqrt(sum((UT - AP(m,:)).^2,2));
        %     Fixed distance
        %     dist = 1 * ones(1,K);
%         Beta = 10.^((-128.1 - 37.6*log10(dist))/10);
        %     Pathloss for normalized noise (variance=1)
%         Beta = Beta/sigma_sqr;
        %     SNR = -10 * log10(0.2) + 25;
        %     Beta = 10^(SNR/10) * ones(1,K);
        %     Calculate gamma_n=a_n*sqrt(beta_n)
%         Beta2 = 4*zeros(K,M);
        gamma(idx,i,m) = Beta2(idx,m);

        % Received signal
        H = (1/sqrt(2)*complex(randn(K,N),randn(K,N)));

        % Rician fading
        N_rician = randperm(K);
        Percent_rician = 0.3;
        K_rician = 0.6*rand(K*Percent_rician,1);
        random_phase_shift = sqrt(-1)*2*pi*rand(K*Percent_rician,1);
        H(N_rician(1:Percent_rician*K),:) = sqrt(K_rician./(K_rician+1)).*exp(random_phase_shift * linspace(0,N-1,N)) + sqrt(1./(K_rician+1)) .* H(N_rician(1:Percent_rician*K),:);
        
        W = 10^(0.1*sqrt(0.2)*randn(1))*(1/sqrt(2)*complex(randn(L,N),randn(L,N)));
%         X = diag(sqrt(gamma(:,i,m)));
        X = diag(sqrt(L*gamma(:,i,m).*10.^(0.1*sqrt(0)*randn(K,1))));
        Y_real(:,:,m) = S*X*H;
        Y(:,:,m) = Y_real(:,:,m) + W;
        % Record G
        G_real(:,:,m,i) = X * H;
        %     Record pathloss5
        pathloss = Beta2';
        %     pathloss_real = Beta_noised';
    end
    gamma_real(:,i) = max(gamma(:,i,:),[],3);
%     maxBeta = max(Beta_real * sigma_sqr);
    fprintf('Set %d\n', i); 
%% Activity detection and Channel estimation
% for m=1:M
    % [~,G_hat_rowAMP(:,:,m,i),~,~,~] = noisyCAMPmmseforKLS(S,K,N,L,Y(:,:,m),G_real(:,:,m,i),100,0.05,Beta'*L,1,1);
%     for j=1:N
        % G_hat_gmmvamp(:,j,m,i) = CAMP(Y(:,j,m),S,1,Beta2(:,m)*L,0.05);
%     end
% end

% GHVI algorithm
% tic;
[G_hat_vigh(:,:,:,i),~] = VIAD_GH_cellfree(Y,S);
% [G_hat_vigh_v2(:,:,:,i)] = VIAD_GH_cellfree_v2(Y,S);
% [G_hat_vigh_v3(:,:,:,i)] = VIAD_GH_cellfree_v3(Y,S);
% runtime_vigh = runtime_vigh + toc/monte;

% [G_hat_vigh(:,:,:,i), z_hat_vigh] = MAP_GH_cellfree(Y,S);
% [G_hat_vigh(:,:,:,i), ~, ~] = VIAD_GG_cellfree(Y,S);
% tic;
% GMMV-AMP
% for m=1:M
% %         [G_hat_vigh(:,:,m,i), ~] = SBLAD(Y(:,:,m),S,M);
    % [G_hat_gmmvamp(:,:,m,i), lambda_gmmvamp(:,:,m,i)] = gmmv_amp(Y(:,:,m), S, 0.3, 50, 1e-5, 0);
    % G_hat_gmmvamp(:,:,m,i) = LSCE(Y(:,:,m),S);
    % [G_hat_CVAMP(:,:,m,i),~,~] = CVAMP_singlecell(Y(:,:,m),S,1,Beta2(:,m)*L,epsilon);
% end
% runtime_map = runtime_map + toc/monte;
% Gaussian-gamma
% tic;
% for m=1:M
%     [G_hat_gmmvamp(:,:,m,i), ~] = VIAD_GG_cellfree(Y(:,:,m),S);
% end
% runtime_vigh = runtime_vigh + toc/monte;

% CAMP
% tic;
% [G_hat_CAMP(:,:,:,i),P_a_CAMP(:,i),conv_CAMP(i)] = CAMP(Y,S,1,Beta2*L,epsilon);
% [G_hat_CVAMP(:,:,:,i),P_a_CVAMP(:,i),conv_CVAMP(i)] = CVAMP(Y,S,1,Beta2*L,epsilon);
% runtime_vigh = runtime_vigh + toc/monte;


end
%% Estimation END
 fprintf('Simulation Finished\n');
%% Best AP
% [N,M,K,monte] = size(G_hat_vigh);
% Ga_hat     = zeros(N,M,monte);
% Ga_real    = zeros(N,M,monte);
% maxbs     = zeros(N,monte);
% 
% for j=1:monte
%     Gnorm = G_hat_vigh(:,:,:,j).*conj(G_hat_vigh(:,:,:,j));
%     normsum = sum(Gnorm,2);
%     [~,maxbs(:,j)] = max(normsum,[],3);
% end
% 
% for j=1:monte
%     for n=1:N
%         Ga_hat(n,:,j) = G_hat_vigh(n,:,maxbs(n,j),j);
%         Ga_real(n,:,j) = G_real(n,:,maxbs(n,j),j); 
%     end
% end
% Ganorm_real = sum(Ga_real .* conj(Ga_real),2);
% Ganorm_hat = sum(Ga_hat .* conj(Ga_hat),2);

%% GMMV-AMP
% G_hat_gmmvamp = G_hat_CAMP;
% G_hat_gmmvamp = G_hat_CVAMP;
G_hat_gmmvamp = G_hat_vigh;
[N,M,K,monte] = size(G_hat_gmmvamp);
Ga_hat_amp     = zeros(N,M,monte);
Ga_real_amp    = zeros(N,M,monte);
maxbs_amp     = zeros(N,monte);

for j=1:monte
    % Gnorm_amp = G_hat_gmmvamp(:,:,:,j).*conj(G_hat_gmmvamp(:,:,:,j));
    Gnorm_amp = G_real(:,:,:,j).*conj(G_real(:,:,:,j));
    normsum_amp = sum(Gnorm_amp,2);
    [~,maxbs_amp(:,j)] = max(normsum_amp,[],3);
end

for j=1:monte
    for n=1:N
        Ga_hat_amp(n,:,j) = G_hat_gmmvamp(n,:,maxbs_amp(n,j),j);
        Ga_real_amp(n,:,j) = G_real(n,:,maxbs_amp(n,j),j); 
    end
end
Ganorm_real_amp = sum(Ga_real_amp .* conj(Ga_real_amp),2);
Ganorm_hat_amp = sum(Ga_hat_amp .* conj(Ga_hat_amp),2);


% %% v2
% Ga_hat_v2     = zeros(N,M,monte);
% Ga_real_v2    = zeros(N,M,monte);
% maxbs_v2     = zeros(N,monte);
% 
% for j=1:monte
%     Gnorm_v2 = G_hat_vigh_v2(:,:,:,j).*conj(G_hat_vigh_v2(:,:,:,j));
%     normsum_v2 = sum(Gnorm_v2,2);
%     [~,maxbs_v2(:,j)] = max(normsum_v2,[],3);
% end
% 
% for j=1:monte
%     for n=1:N
%         Ga_hat_v2(n,:,j) = G_hat_vigh_v2(n,:,maxbs(n,j),j);
%         Ga_real_v2(n,:,j) = G_real(n,:,maxbs(n,j),j); 
%     end
% end
% Ganorm_real_v2 = sum(Ga_real_v2 .* conj(Ga_real_v2),2);
% Ganorm_hat_v2 = sum(Ga_hat_v2 .* conj(Ga_hat_v2),2);
% 
% %% v3
% Ga_hat_v3     = zeros(N,M,monte);
% Ga_real_v3    = zeros(N,M,monte);
% maxbs_v3     = zeros(N,monte);
% 
% for j=1:monte
%     Gnorm_v3 = G_hat_vigh_v3(:,:,:,j).*conj(G_hat_vigh_v3(:,:,:,j));
%     normsum_v3 = sum(Gnorm_v3,2);
%     [~,maxbs_v3(:,j)] = max(normsum_v3,[],3);
% end
% 
% for j=1:monte
%     for n=1:N
%         Ga_hat_v3(n,:,j) = G_hat_vigh_v3(n,:,maxbs(n,j),j);
%         Ga_real_v3(n,:,j) = G_real(n,:,maxbs_v3(n,j),j); 
%     end
% end
% Ganorm_real_v3 = sum(Ga_real_v3 .* conj(Ga_real_v3),2);
% Ganorm_hat_v3 = sum(Ga_hat_v3 .* conj(Ga_hat_v3),2);