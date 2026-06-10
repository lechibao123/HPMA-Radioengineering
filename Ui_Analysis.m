function F = Ui_Analysis(x, Ns, m1a, xi1)
% CDF_TAS1_FUV: Calculates the Cumulative Distribution Function (CDF) for the FUV
%               channel according to equation (26) from the paper "Performance
%               Analysis of NOMA-Enabled Vehicular Communication Systems...".
%
% Syntax:
%   F = CDF_TAS1_FUV(x, Ns, m1a, xi1)
%
% INPUTS:
%   x   - Channel gain value. Can be a scalar or a vector.
%   Ns  - Number of transmit antennas (positive integer).
%   m1a - Fading parameter 'm' for the S -> FUV link (real number, typically >= 0.5).
%   xi1 - Channel parameter xi_1 (positive real number).
%
% OUTPUT:
%   F   - The value of the Cumulative Distribution Function (CDF) at points x.
%
% Reference:
%   [1] Jaiswal, N., & Purohit, N. (2021). Performance Analysis of
%       NOMA-Enabled Vehicular Communication Systems With Transmit Antenna
%       Selection Over Double Nakagami-m Fading. IEEE Transactions on

% --- 1. Input Validation ---
if any(x < 0)
    error('Input value x (channel gain) cannot be negative.');
end
if Ns <= 0 || mod(Ns, 1) ~= 0
    error('Number of antennas Ns must be a positive integer.');
end
if m1a < 0.5
    warning('A value of m1a < 0.5 might lead to unexpected results.');
end

% --- 2. Calculate Intermediate Parameters ---
% According to the paper, |m1b - m1a| = 1/2. [cite_start]We choose m1b = m1a + 0.5[cite: 296].
m1b = m1a + 0.5;

% Calculate chi_1 according to the definition below eq. (206) [cite_start][cite: 206].
chi1 = (2^(1-2*m1a) * sqrt(pi)) / (gamma(m1a) * gamma(m1b));

multiplier = (chi1 * gamma(2*m1a))^Ns;

% --- 3. Calculate p_r and q_r^n Coefficients ---
p_limit = floor(2*m1a - 1);
p_coeffs = zeros(1, p_limit + 1);
for r_idx = 0:p_limit
    p_coeffs(r_idx + 1) = (1/factorial(r_idx)) * (4*xi1)^(r_idx/2);
end

% Calculate q_r^n coefficients using dynamic programming to avoid recalculation.
q_max_r = floor(Ns * (2*m1a - 1));
Q = zeros(Ns + 1, q_max_r + 1);

for n = 0:Ns
    Q(n + 1, 1) = p_coeffs(1)^n;
    
    n_limit_r = floor(n * (2*m1a - 1));
    if n > 0 && n_limit_r > 0
        for r = 1:n_limit_r
            sum_val = 0;
            limit_l = min(r, p_limit);
            for l = 1:limit_l
                sum_val = sum_val + (l*n - r + l) * p_coeffs(l+1) * Q(n+1, r-l+1);
            end
            % p_coeffs(1) is p_0, and p_0 = 1, so we divide by r.
            Q(n+1, r+1) = (1/r) * sum_val;
        end
    end
end

% --- 4. Calculate the Main Summation of Equation (26) ---
% Ensure x is a column vector for efficient matrix operations.
x_col = x(:);
total_sum = zeros(size(x_col));

for n = 0:Ns
    n_limit_r = floor(n * (2*m1a - 1));
    for r = 0:n_limit_r
        q_val = Q(n + 1, r + 1);
        
        term = nchoosek(Ns, n) * (-1)^n * q_val ...
             .* (x_col.^(r/2)) ...
             .* exp(-2*n*sqrt(x_col*xi1));
        
        total_sum = total_sum + term;
    end
end

% --- 5. Final Result ---
F = multiplier * total_sum;

% Reshape the output to match the original size of input x.
F = reshape(F, size(x));

end