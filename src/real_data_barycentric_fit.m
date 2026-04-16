clc;clear; 
data = load('20260224_0+Y,1Y,2Y,3-Y,4-Y,5Y,6XYZ.mat');
% %% Load data file
% [fileName, folderPath] = uigetfile('*.mat', 'Select data file');
% if fileName ~= 0
%     fullPath = fullfile(folderPath, fileName);
%     data = load(fullPath);
% else
%     disp('No file selected');
%     return;
% end
%acc = zeros(size(data.mtxEU(:, 2, 1)));
acc = data.mtxEU(:,2, 1);
%Voltage = data.mtxVoltage(:, 3)

%for i = 1:3
%    temp = data.mtxEU(:, 2, i);
%    acc = acc + temp;
%    plot(acc);
%end
%acc = acc/3;
% User inputs
Fs = 1707;          % sampling frequency in Hz (set to your value)
acc = acc(:);       % ensure column vector

% Design parameters
fc = 70;            % cutoff frequency in Hz 從0濾到120
n = 4;              % filter order

% Normalized cutoff (0..1 where 1 = Nyquist = Fs/2)
Wn = fc / (Fs/2);

% Option A: Simple design (b,a) and zero-phase filtering
[b,a] = butter(n, Wn, 'low');        % design IIR
acc_filt_A = filtfilt(b, a, acc);    % zero-phase forward-backward
fs = 1707;
% FFT (complex spectrum)
N   = numel(acc);
Y   = fft(acc);
f   = (0:floor(N/2)).' * (fs/N);           % single-sided freq axis
Hc  = Y(1:floor(N/2)+1) / N;               % complex spectrum (scaled)

% pick frequency range to fit (like your example)
f_min = 3;
f_max = 9;
idx   = (f >= f_min) & (f <= f_max);

f_fit = f(idx);
Hdata = Hc(idx);                           % <<< complex data (recommended)
% If you only want magnitude fitting, you can do: Hdata = abs(Hc(idx));

fprintf('Data generated: %d frequency points to fit\n', numel(f_fit));

% Laplace variable on imaginary axis
s_all = 1i*2*pi*f_fit;

%% ============================================================
% PART 2) Greedy support point selection + barycentric fit
% ============================================================

Nsupport_target   = 9 ;      % Fig.13 style
Ninit             = 3;       % initial supports
exclude_radius_hz = 0.001;     % avoid selecting too close (Hz)
lambda_reg        = 1e-8;    % small regularization (stabilizes SVD)

% noise floor (Hugh: ~1/1000 of data)
noise_floor_abs = max(abs(Hdata))/1000;
noise_floor_rel = 1/1000;

% init supports: edges + peak (more "Hugh-like" than even spacing)
[~,iL] = min(f_fit);
[~,iR] = max(f_fit);
[~,iP] = max(abs(Hdata));
idx_support = unique([iL; iR; iP]);

% greedy iterations
Hfit = zeros(size(Hdata));
relErr = zeros(size(Hdata));

maxIters = max(0, Nsupport_target - numel(idx_support));
for it = 1:maxIters

    sn = s_all(idx_support);
    hn = Hdata(idx_support);

    % (A) linear weight solve via SVD null-space
    w = solve_weights_svd(s_all, Hdata, sn, hn, lambda_reg);

    % (B) evaluate barycentric model
    Hfit = eval_barycentric(s_all, sn, hn, w);

    % (C) relative error (stable denominator)
    absErr = abs(Hfit - Hdata);
    relErr = absErr ./ max(abs(Hdata), noise_floor_abs);

    % (D) greedy pick: max relative error, excluding nearby supports
    idx_new = pick_next_support(f_fit, absErr, idx_support, exclude_radius_hz);
    idx_support(end+1,1) = idx_new; %#ok<SAGROW>

    fprintf('Iter %2d | #support=%2d | max(relErr)=%.2e at %.2f Hz\n', ...
        it, numel(idx_support), max(absErr), f_fit(idx_new));

    % optional early stop
    if max(relErr) < noise_floor_rel
        fprintf('Stop early: reached noise floor.\n');
        break;
    end
end

% final evaluation
sn = s_all(idx_support);
hn = Hdata(idx_support);
w  = solve_weights_svd(s_all, Hdata, sn, hn, lambda_reg);
Hfit  = eval_barycentric(s_all, sn, hn, w);
absErr = abs(Hfit - Hdata);
relErr = absErr ./ max(abs(Hdata), noise_floor_abs);

%% ============================================================
% PART 3) Plot (single 2D figure like Fig.13)
% ============================================================

figure('Color','w'); hold on; grid on;

yyaxis left
set(gca,'YScale','log')
semilogy(f_fit, abs(Hdata), 'b', 'LineWidth',2);
semilogy(f_fit, abs(Hfit),  'r--','LineWidth',2);
semilogy(f_fit(idx_support), abs(Hdata(idx_support)), 'ko', ...
    'MarkerFaceColor','k','MarkerSize',6);
ylabel('|H|',FontSize=14);

yyaxis right
set(gca,'YScale','log', fontsize= 12)
semilogy(f_fit, absErr, 'g', 'LineWidth',1.2);
yline(noise_floor_rel,'--m','1/1000 noise floor','LineWidth',1.2);
ylabel('Absolute error  |Hfit-Hdata| / |Hdata|');

xlabel('Frequency (Hz)',FontSize=14);
title('Barycentric rational fit (Duffing ODE45 data) + greedy support points', FontSize=18);
legend('Data','Fit','Support points','Abs. error','Noise floor', ...
       'Location','southwest',fontsize = 16);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function w = solve_weights_svd(s_all, Hdata, sn, hn, lambda)
% Build A(k,n) = (hn(n) - Hdata(k)) / (s_k - sn(n))
% Want A*w ≈ 0 -> take last right singular vector of A

    K = numel(s_all);
    N = numel(sn);

    A = zeros(K,N);
    for n = 1:N
        A(:,n) = (hn(n) - Hdata) ./ (s_all - sn(n));
    end

    % remove rows at exact support points (avoid Inf)
    mask = true(K,1);
    for n = 1:N
        mask = mask & (abs(s_all - sn(n)) > 1e-12);
    end
    A = A(mask,:);

    % regularization
    if lambda > 0
        A = [A; sqrt(lambda)*eye(N)];
    end

    [~,~,V] = svd(A,'econ');
    w = V(:,end);

    % normalize (scale cancels in barycentric formula)
    w = w ./ max(abs(w));
end

function H = eval_barycentric(s_all, sn, hn, w)
% H(s) = (sum w_n*h_n/(s-sn)) / (sum w_n/(s-sn))

    K = numel(s_all);
    N = numel(sn);

    numer = zeros(K,1);
    denom = zeros(K,1);

    for n = 1:N
        diff  = s_all - sn(n);
        numer = numer + (w(n)*hn(n))./diff;
        denom = denom +  w(n)./diff;
    end

    H = numer ./ denom;

    % exact hits: enforce interpolation
    for n = 1:N
        hit = abs(s_all - sn(n)) < 1e-12;
        H(hit) = hn(n);
    end
end

function idx_new = pick_next_support(f, metric, idx_support, radius_hz)
% pick argmax(metric) excluding points too close to existing supports

    score = metric(:);
    valid = true(size(f));

    for k = 1:numel(idx_support)
        valid = valid & (abs(f - f(idx_support(k))) > radius_hz);
    end
    score(~valid) = -Inf;

    [~, idx_new] = max(score);

    % fallback
    if ~isfinite(score(idx_new))
        [~, idx_new] = max(metric);
    end
end
