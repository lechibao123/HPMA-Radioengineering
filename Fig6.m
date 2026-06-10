clear; clc; tic;

%% Parameters %%
% Bandwidth allocation
alpha_1 = 0.45; alpha_2 = 0.45; alpha_12 = 0.1;

% Power allocation
delta_HPMA_1 = 0.4; delta_HPMA_2 = 1 - delta_HPMA_1;
delta_NOMA_2 = 0.85; delta_NOMA_1 = 1 - delta_NOMA_2;

% Number Antennas
K = 1;

% Target rates bps/Hz

SNR_dB = 0:0.15:3.9;
R1_arr = SNR_dB;


% Nakagami-m parameters
m_1a = 2;
m_2a = m_1a;

m_1b = m_1a + 0.5;
m_2b = m_2a + 0.5;


% Average power per fading link
Omega_1 = 1;
Omega_2 = 1;

Omega_1a = Omega_1;
Omega_1b = Omega_1;
Omega_2a = Omega_2;
Omega_2b = Omega_2;

%_
zeta = 2;
L = 1;
H = 1;
r = 0.1;
O = 0;
d1 = sqrt(H^2 + r^2 + L^2 - 2*r*L*cosd(O));
d2 = sqrt(H^2 + r^2 + L^2 + 2*r*L*cosd(O));

%_

rho = db2pow(12);

%_
beta_1 = 1 / (alpha_1 + alpha_12);
beta_2 = 1 / (alpha_2 + alpha_12);

% Init
leng_colum = length(R1_arr);
P_out_HPMA_1_sim = zeros(1,leng_colum); P_out_HPMA_2_sim = zeros(1,leng_colum);
P_out_NOMA_1_sim = zeros(1,leng_colum); P_out_NOMA_2_sim = zeros(1,leng_colum);
P_out_HPMA_1 = zeros(1,leng_colum); P_out_HPMA_2 = zeros(1,leng_colum);
P_out_NOMA_1 = zeros(1,leng_colum); P_out_NOMA_2 = zeros(1,leng_colum);
P_out_HPMA_1_asym = zeros(1,leng_colum); P_out_HPMA_2_asym = zeros(1,leng_colum);
P_out_NOMA_1_asym = zeros(1,leng_colum); P_out_NOMA_2_asym = zeros(1,leng_colum); 
LOOP = 1e7;

% Generate random Nakagami-m fading
h_1i_a = gamrnd(m_1a, Omega_1a/m_1a, K, LOOP);
h_1i_b = gamrnd(m_1b, Omega_1b/m_1b, K, LOOP);
h_2i_a = gamrnd(m_2a, Omega_2a/m_2a, K, LOOP);
h_2i_b = gamrnd(m_2b, Omega_2b/m_2b, K, LOOP);

eta = 1;

H_1i = eta * d1^-zeta .* max(h_1i_a .* h_1i_b,[],1);
H_2i = eta * d2^-zeta .* max(h_2i_a .* h_2i_b,[],1);


%% Loop over SNR
for idx = 1:length(R1_arr)
    R1 = R1_arr(idx);
    R2 = R1;
    gamma_th_1 = 2^R1 - 1;
    gamma_th_2 = 2^R2 - 1;

    
    %% HPMA SINRs
    gamma_HPMA_1 = ((alpha_1 + alpha_12*delta_HPMA_1) * beta_1 * rho .* H_1i) ./ ...
        (alpha_12 * delta_HPMA_2 * beta_1 * rho .* H_1i + 1);
    gamma_HPMA_2 = ((alpha_2 + alpha_12*delta_HPMA_2) * beta_2 * rho .* H_2i) ./ ...
        (alpha_12 * delta_HPMA_1 * beta_2 * rho .* H_2i + 1);
    
    %% NOMA SINRs
    gamma_NOMA_2 = (delta_NOMA_2 * rho .* H_2i) ./ ...
        (delta_NOMA_1 * rho .* H_2i + 1);
    gamma_NOMA_1_x2 = (delta_NOMA_2 * rho .* H_1i) ./ ...
        (delta_NOMA_1 * rho .* H_1i + 1);
    gamma_NOMA_1_x1 = delta_NOMA_1 * rho .* H_1i;
    
    %% Simulated Outage
    P_out_HPMA_1_sim(idx) = (1 -  mean(gamma_HPMA_1 < gamma_th_1)).*R1;
    P_out_HPMA_2_sim(idx) = (1 - mean(gamma_HPMA_2 < gamma_th_2)).*R2;
    P_out_NOMA_1_sim(idx) = (1 -  mean((gamma_NOMA_1_x2 < gamma_th_2) | (gamma_NOMA_1_x1 < gamma_th_1))).*R1;
    P_out_NOMA_2_sim(idx) = (1 - mean(gamma_NOMA_2 < gamma_th_2)).*R2;
    

    %% Analytical HPMA
    Lambda1 = gamma_th_1 / (rho * eta * d1^-zeta * beta_1 * ((alpha_1 + alpha_12*delta_HPMA_1) - gamma_th_1 * alpha_12 * delta_HPMA_2));
    Lambda2 = gamma_th_2 / (rho * eta * d2^-zeta * beta_2 * ((alpha_2 + alpha_12*delta_HPMA_2) - gamma_th_2 * alpha_12 * delta_HPMA_1));
    
    %% Analytical NOMA
    Lambda3 = gamma_th_2 / (rho * eta * d2^-zeta * (delta_NOMA_2 - gamma_th_2 * delta_NOMA_1));
    
    
    Lambda4 = gamma_th_2 / (rho * eta * d1^-zeta * (delta_NOMA_2 - gamma_th_2 * delta_NOMA_1));
    Lambda5 = gamma_th_1 / (rho * eta * d1^-zeta * delta_NOMA_1);
    Lambda_max = max(Lambda4,Lambda5);
    
    xi_1 = (m_1a*m_1b)/Omega_1a/Omega_1b;
    xi_2 = (m_2a*m_2b)/Omega_2a/Omega_2b;
    
    chi_1 = sqrt(pi) * 2^(1 - 2*m_1a) * (1/gamma(m_1a)/gamma(m_1b));
    chi_2 = sqrt(pi) * 2^(1 - 2*m_2a) * (1/gamma(m_2a)/gamma(m_2b));

    if Lambda1 < 0
        P_out_HPMA_1(idx) = P_out_HPMA_1_sim(idx);
    else
        P_out_HPMA_1(idx) = (1 - Ui_Analysis(Lambda1, K, m_1a, xi_1)).*R1;
    end
    
    if Lambda2 < 0
       P_out_HPMA_2(idx) = P_out_HPMA_2_sim(idx);
    else
       P_out_HPMA_2(idx) = (1 - Ui_Analysis(Lambda2, K, m_2a, xi_2)).*R2;
    end
        
    
    
    P_out_NOMA_1(idx) = (1 - Ui_Analysis(Lambda_max, K, m_1a, xi_1)).*R1;
%     P_out_NOMA_2(idx) = (1 - Ui_Analysis(Lambda3, K, m_2a, xi_2)).*R2;
end


%% Plot %%
blue1 = [0.00,0.45,0.74];  pink1 = [1.00,0.07,0.65];
green1 = [0.47,0.67,0.19]; orrange = [0.85,0.33,0.10];

plot(SNR_dB,P_out_HPMA_1,':o','Linewidth',1,'MarkerFaceColor','w','MarkerSize',9,'Color',green1); hold on
plot(SNR_dB,P_out_HPMA_2,':s','Linewidth',1,'MarkerFaceColor','w','MarkerSize',9,'Color',pink1);
plot(SNR_dB,P_out_NOMA_1_sim,':d','Linewidth',1,'MarkerFaceColor','w','MarkerSize',9,'Color',blue1);
plot(SNR_dB,P_out_NOMA_2_sim,':v','Linewidth',1,'MarkerFaceColor','w','MarkerSize',9,'Color','m');

plot(SNR_dB,P_out_HPMA_1_sim,'k.','Linewidth',1,'MarkerSize',9);
plot(SNR_dB,P_out_HPMA_2_sim,'k.','Linewidth',1,'MarkerSize',9);
plot(SNR_dB,P_out_NOMA_1_sim,'k.','Linewidth',1,'MarkerSize',9);
plot(SNR_dB,P_out_NOMA_2_sim,'k.','Linewidth',1,'MarkerSize',9);


xlabel('$R_1 = R_2$ [bps/Hz]','Interpreter','latex','fontsize',15);
ylabel('Throughput','Interpreter','latex','fontsize',15);
h = legend({'User 1, HPMA','User 2, HPMA','User 1, NOMA','User 2, NOMA','Simulation'},'Interpreter','latex','Location','best');
set(h,'interpreter','latex','fontsize',11);

ylim([0 1.7]);
toc;