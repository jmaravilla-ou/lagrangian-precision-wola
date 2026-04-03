% e_listen - informal listening test, not part of the paper results pipeline.
%
% picks one random clip from dev-clean and one from dev-other, runs each
% through the full wola signal path under all three allocation strategies
% at all four budget levels, and writes labelled wav files so you can
% compare them by ear. also prints snr, stoi, D_w, and D_mse per condition.
%
% signal path: peak-normalize -> hann window -> 128-pt fft ->
%   operand isolation mask -> midrise quantizer -> 128-pt ifft ->
%   synthesis window -> overlap-add
%
% outputs (written to ./E_listen_output/):
%   <clip>_reference.wav
%   <clip>_uniform_<pct>pct.wav
%   <clip>_linear_<pct>pct.wav
%   <clip>_proposed_<pct>pct.wav
%   E_listen_report.csv
%
% requires: matlab audio toolbox (stoi, audiowrite)
%           allocators.m on path
%           E1_sigma2_raw.csv and E1_sigma2_eff.csv in working directory
clear; clc; close all;

%% ── Parameters ──────────────────────────────────────────────────────────
LIBRISPEECH_ROOTS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'
};

M           = 64;
FS_TARGET   = 16000;
B_MAX       = 24;
B_MIN       = 4;
ALPHA       = 1.37;
BETA        = 0.00;
FRAME       = 2*M;          % 128 samples
HOP         = M;            % 64 samples — 50% overlap
WIN         = hann(FRAME, 'periodic');
C_TOTAL_REF = M * 16;       % 1024 bits
BUDGET_FRACS = [0.40, 0.60, 0.80, 1.00];
BUDGET_PCTS  = [40,   60,   80,   100  ];

% Output folder
OUT_DIR = 'E_listen_output';
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ── Psychoacoustic weight (must match every other experiment) ────────────
fc = (1:M)' * (FS_TARGET/2) / M;   % 125, 250, ..., 8000 Hz
S  = exp(-0.5 * (log(fc / 2500) / 0.55).^2);
S  = S / max(S);

%% ── Load sigma2 vectors from E1 ─────────────────────────────────────────
% sigma2_raw: used for D_w and D_mse evaluation
% sigma2_eff: used as allocator input
if ~isfile('E1_sigma2_raw.csv') || ~isfile('E1_sigma2_eff.csv')
    error(['E1_sigma2_raw.csv or E1_sigma2_eff.csv not found.\n' ...
           'Run E1 first to generate these files.']);
end
Tr = readtable('E1_sigma2_raw.csv', 'VariableNamingRule', 'preserve');
Te = readtable('E1_sigma2_eff.csv', 'VariableNamingRule', 'preserve');
sigma2_raw = Tr.sigma2_raw;
sigma2_eff = Te.sigma2_eff;
fprintf('Loaded sigma2_raw and sigma2_eff (%d bands).\n\n', M);

%% ── Discover corpus and pick one random clip per corpus ─────────────────
corpus = discover_librispeech(LIBRISPEECH_ROOTS);
N = numel(corpus);
if N == 0, error('No .flac files found. Check LIBRISPEECH_ROOTS.'); end

% Separate dev-clean and dev-other
is_clean = strcmp({corpus.corpus}, 'dev-clean');
is_other = strcmp({corpus.corpus}, 'dev-other');

rng('shuffle');   % different clip every run
picks = {};
if any(is_clean)
    idx = find(is_clean);
    picks{end+1} = corpus(idx(randi(numel(idx))));
end
if any(is_other)
    idx = find(is_other);
    picks{end+1} = corpus(idx(randi(numel(idx))));
end
if isempty(picks)
    picks{1} = corpus(randi(N));
end

fprintf('Selected clips:\n');
for p = 1:numel(picks)
    fprintf('  [%d] %s  (%s)\n', p, picks{p}.label, picks{p}.corpus);
end
fprintf('\n');

%% ── Pre-compute allocations (same for all clips) ─────────────────────────
strategies   = {'Uniform', 'Linear', 'Proposed'};
N_STRAT      = 3;
N_BUD        = numel(BUDGET_FRACS);
B_alloc      = cell(N_BUD, N_STRAT);

fprintf('Computing allocations...\n');
for b = 1:N_BUD
    C_total      = round(C_TOTAL_REF * BUDGET_FRACS(b));
    B_alloc{b,1} = allocate_uniform(C_total, M);
    B_alloc{b,2} = allocate_linear(sigma2_eff, C_total, M, B_MIN, B_MAX);
    B_alloc{b,3} = allocate_quadratic(sigma2_eff, ALPHA, BETA, C_total, M, B_MIN, B_MAX);
    fprintf('  %3d%%: Uniform sum=%d  Linear sum=%d  Proposed sum=%d\n', ...
            BUDGET_PCTS(b), sum(B_alloc{b,1}), sum(B_alloc{b,2}), sum(B_alloc{b,3}));
end
fprintf('\n');

%% ── Print per-band allocation at 40% for inspection ─────────────────────
fprintf('Per-band allocation at 40%% budget:\n');
fprintf('%-6s  %-8s  %-4s  %-4s  %-4s  %-4s\n', ...
        'Band', 'Freq(Hz)', 'S(k)', 'Unif', 'Lin', 'Prop');
fprintf('%s\n', repmat('-',1,44));
for k = 1:M
    if mod(k-1,4)==0 || k==1   % print every 4th band + all first few
        fprintf('  k=%-3d  %-8.0f  %.3f  %-4d  %-4d  %-4d\n', ...
                k, fc(k), S(k), ...
                B_alloc{1,1}(k), B_alloc{1,2}(k), B_alloc{1,3}(k));
    end
end
fprintf('\n');

%% ── Main loop: process each selected clip ────────────────────────────────
report = {};   % collect rows for CSV

for p = 1:numel(picks)
    clip = picks{p};
    fprintf('========================================\n');
    fprintf('Clip: %s  (%s)\n', clip.label, clip.corpus);
    fprintf('========================================\n');

    % Load and prepare
    [x_raw, fs] = audioread(clip.path);
    x_raw = mean(x_raw, 2);                        % mono
    if fs ~= FS_TARGET
        x_raw = resample(x_raw, FS_TARGET, fs);
    end
    x_ref = x_raw / (max(abs(x_raw)) + eps);       % peak-normalise to [-1,1]
    sig_len = length(x_ref);

    % Write reference
    ref_name = fullfile(OUT_DIR, sprintf('%s_reference.wav', clip.label));
    audiowrite(ref_name, x_ref, FS_TARGET);
    fprintf('  Reference written: %s\n', ref_name);

    % WOLA analysis (complex — needed for resynthesis)
    [Xk_mag, Xk_phase] = wola_analyze_complex(x_ref, M, WIN, HOP);

    fprintf('\n  %-12s  %-6s  %-10s  %-10s  %-10s  %-10s  %-8s\n', ...
            'Strategy', 'Budget', 'SNR(dB)', 'STOI', 'D_w', 'D_mse', 'D_w/D_mse');
    fprintf('  %s\n', repmat('-', 1, 72));

    for b = 1:N_BUD
        pct = BUDGET_PCTS(b);
        for s = 1:N_STRAT
            Bv = B_alloc{b,s};

            % ── Operand isolation + quantization ─────────────────────────
            % Models the circuit path: AND-mask zeros the (B_max - Bk) LSBs
            % of both operands before the multiply. For the magnitude subband
            % representation, this is equivalent to a midrise uniform quantizer
            % with 2^(Bk-1) steps over the signal's dynamic range.
            Xk_q = operand_iso_quantize(Xk_mag, Bv, B_MAX);

            % ── Resynthesis: IFFT + synthesis window + OLA ───────────────
            y = wola_synthesize(Xk_q, Xk_phase, M, WIN, HOP, sig_len);
            y = max(-1, min(1, y));                 % clip guard

            % ── Metrics ──────────────────────────────────────────────────
            snr_val  = compute_snr(x_ref, y);
            stoi_val = stoi(x_ref, y, FS_TARGET);

            % Analytical D_w and D_mse using the actual per-band allocation
            dw_val   = sum(S        .* sigma2_raw .* 2.^(-2*Bv));
            dmse_val = sum(sigma2_raw              .* 2.^(-2*Bv));
            ratio    = dw_val / (dmse_val + eps);

            % Empirical per-clip per-band MSE (what actually happened)
            band_mse = compute_band_mse(Xk_mag, Xk_q);
            dw_emp   = sum(S .* band_mse);
            dm_emp   = sum(band_mse);

            strat_name = strategies{s};
            fprintf('  %-12s  %3d%%   %8.2f    %.4f   %.4e  %.4e  %.2e\n', ...
                    strat_name, pct, snr_val, stoi_val, dw_val, dmse_val, ratio);

            % Write audio
            fname = fullfile(OUT_DIR, sprintf('%s_%s_%dpct.wav', ...
                             clip.label, lower(strat_name), pct));
            audiowrite(fname, y, FS_TARGET);

            % Collect report row
            report{end+1} = {clip.label, clip.corpus, strat_name, pct, ...
                             snr_val, stoi_val, dw_val, dmse_val, ratio, ...
                             dw_emp, dm_emp, ...
                             sum(Bv), min(Bv), max(Bv)};
        end
        fprintf('\n');
    end

    % Quick ratio commentary
    fprintf('  --- D_w/D_mse ratio interpretation ---\n');
    fprintf('  A ratio << 1 means the perceptual metric heavily discounts\n');
    fprintf('  noise in low-sensitivity bands (low-freq large-sigma bands).\n');
    fprintf('  The absolute D_mse value looks large because low-freq bands\n');
    fprintf('  have high raw energy (sigma2_raw up to 1.0 at 375 Hz) but\n');
    fprintf('  S(k) near 0, so they receive few bits and emit large raw noise.\n');
    fprintf('  If the audio sounds clean, D_mse is misleading; D_w is honest.\n\n');
end

%% ── Write report CSV ─────────────────────────────────────────────────────
fid = fopen(fullfile(OUT_DIR, 'E_listen_report.csv'), 'w');
fprintf(fid, ['Clip,Corpus,Strategy,BudgetPct,SNR_dB,STOI,' ...
              'Dw_analytical,Dmse_analytical,Dw_Dmse_ratio,' ...
              'Dw_empirical,Dmse_empirical,' ...
              'sum_Bk,min_Bk,max_Bk\n']);
for i = 1:numel(report)
    r = report{i};
    fprintf(fid, '%s,%s,%s,%d,%.4f,%.6f,%.6e,%.6e,%.4e,%.6e,%.6e,%d,%d,%d\n', ...
            r{1},r{2},r{3},r{4},r{5},r{6},r{7},r{8},r{9},r{10},r{11},r{12},r{13},r{14});
end
fclose(fid);
fprintf('Report saved: %s\n', fullfile(OUT_DIR,'E_listen_report.csv'));
fprintf('\nAudio files written to: %s/\n', OUT_DIR);
fprintf('Open them in MATLAB, Audacity, or any audio player to listen.\n');
fprintf('Compare _uniform_ vs _linear_ vs _proposed_ at the same budget.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function [Xk_mag, phase] = wola_analyze_complex(x, M, win, hop)
    % 128-point FFT, keep first M=64 bins (positive frequencies).
    frame_len = 2*M;
    n_frames  = floor((length(x) - frame_len) / hop) + 1;
    Xk_mag    = zeros(M, n_frames);
    phase     = zeros(M, n_frames);
    for f = 1:n_frames
        seg       = x((f-1)*hop + 1 : (f-1)*hop + frame_len);
        spec      = fft(seg .* win);
        spec_half = spec(1:M);
        Xk_mag(:,f) = abs(spec_half);
        phase(:,f)  = angle(spec_half);
    end
end


function Xk_q = operand_iso_quantize(Xk_mag, Bv, B_max)
    % Models AND-masking the (B_max - Bk) LSBs of both MAC operands.
    % Equivalent to a midrise uniform quantizer with step 2^(-(Bk-1)).
    % For signed two's-complement: zeroing LSBs preserves the sign bit;
    % the quantization error is bounded by 2^(-(Bk-1)) per sample.
    [M, n_frames] = size(Xk_mag);
    Xk_q = zeros(M, n_frames);
    for k = 1:M
        bk         = max(1, min(B_max, round(Bv(k))));
        scale      = 2^(bk - 1);          % number of quantisation levels / 2
        Xk_q(k,:)  = round(Xk_mag(k,:) * scale) / scale;
    end
end


function x_out = wola_synthesize(Xk_mag, phase, M, win, hop, sig_len)
    % Reconstruct real signal via conjugate-symmetric IFFT + overlap-add.
    frame_len = 2*M;
    n_frames  = size(Xk_mag, 2);
    x_out     = zeros(sig_len + frame_len, 1);
    norm_acc  = zeros(sig_len + frame_len, 1);
    for f = 1:n_frames
        spec_half = Xk_mag(:,f) .* exp(1j * phase(:,f));
        spec_full = zeros(frame_len, 1);
        spec_full(1)         = real(spec_half(1));   % DC — real
        spec_full(2:M)       = spec_half(2:M);
        spec_full(M+1)       = 0;                    % Nyquist — zero
        spec_full(M+2:end)   = conj(flipud(spec_half(2:M)));
        frame = real(ifft(spec_full));
        frame_win = frame .* win;                    % synthesis window
        idx  = (f-1)*hop + (1:frame_len);
        idx  = idx(idx <= length(x_out));
        x_out(idx)    = x_out(idx)    + frame_win(1:length(idx));
        norm_acc(idx) = norm_acc(idx) + win(1:length(idx)).^2;
    end
    norm_acc(norm_acc < 1e-8) = 1;
    x_out = x_out(1:sig_len) ./ norm_acc(1:sig_len);
end


function snr = compute_snr(ref, deg)
    % SNR in dB: 10*log10(sum(ref^2) / sum((ref-deg)^2))
    noise = ref - deg;
    snr   = 10 * log10(sum(ref.^2) / (sum(noise.^2) + eps));
end


function band_mse = compute_band_mse(Xk_ref, Xk_q)
    % Per-band mean-squared error between reference and quantised magnitudes.
    M        = size(Xk_ref, 1);
    band_mse = zeros(M, 1);
    for k = 1:M
        band_mse(k) = mean((Xk_ref(k,:) - Xk_q(k,:)).^2);
    end
end


function Bv = allocate_uniform(C_total, M)
    B_base   = floor(C_total / M);
    rem_bits = C_total - B_base * M;
    Bv       = B_base * ones(M, 1);
    Bv(1:rem_bits) = Bv(1:rem_bits) + 1;
end

function Bv = allocate_linear(sigma2, C_total, M, B_floor, B_ceil)
    Bv = zeros(M, 1);
    lo = 1e-20; hi = max(sigma2) * 1e8;
    for iter = 1:600
        lam  = (lo + hi) / 2;
        raw  = 0.5 * log2(max(eps, sigma2 / lam));
        Bv   = min(B_ceil, max(B_floor, raw));
        cost = sum(Bv);
        if cost > C_total, lo = lam; else, hi = lam; end
        if abs(cost - C_total) / C_total < 1e-11, break; end
    end
    Bv_int   = min(B_ceil, max(B_floor, round(Bv)));
    leftover = C_total - sum(Bv_int);
    [~, ord] = sort(sigma2, 'descend');
    i = 1;
    while leftover ~= 0 && i <= M
        k = ord(i);
        if leftover > 0 && Bv_int(k) < B_ceil
            Bv_int(k) = Bv_int(k) + 1; leftover = leftover - 1;
        elseif leftover < 0 && Bv_int(k) > B_floor
            Bv_int(k) = Bv_int(k) - 1; leftover = leftover + 1;
        end
        i = i + 1;
    end
    Bv = double(Bv_int);
end

function Bv = allocate_quadratic(sigma2, alpha, beta, C_bits, M, B_floor, B_ceil)
    if nargin < 6, B_floor = 4;  end
    if nargin < 7, B_ceil  = 24; end
    B_avg = C_bits / M;
    C_lut = M * (alpha * B_avg^2 + beta * B_avg);
    Bv = B_floor * ones(M, 1);
    lo = 1e-30; hi = 1e15;
    for iter = 1:600
        lam = (lo + hi) / 2;
        for k = 1:M
            Bv(k) = solve_Bk_newton(sigma2(k), alpha, beta, lam, B_floor);
        end
        cost = sum(alpha * Bv.^2 + beta * Bv);
        if cost > C_lut, lo = lam; else, hi = lam; end
        if abs(cost - C_lut) / C_lut < 1e-11, break; end
    end
    Bv = min(B_ceil, max(B_floor, Bv));
    Bv_int   = min(B_ceil, max(B_floor, round(Bv)));
    leftover = C_bits - sum(Bv_int);
    [~, ord] = sort(sigma2, 'descend');
    i = 1;
    while leftover ~= 0 && i <= M
        k = ord(i);
        if leftover > 0 && Bv_int(k) < B_ceil
            Bv_int(k) = Bv_int(k) + 1; leftover = leftover - 1;
        elseif leftover < 0 && Bv_int(k) > B_floor
            Bv_int(k) = Bv_int(k) - 1; leftover = leftover + 1;
        end
        i = i + 1;
    end
    Bv = double(Bv_int);
end

function Bk = solve_Bk_newton(s2, alpha, beta, lam, B_floor)
    if nargin < 5, B_floor = 4; end
    if s2 <= 0, Bk = B_floor; return; end
    Bk = max(B_floor, 0.5 * log2(max(eps, s2 / (lam + eps))));
    for n = 1:100
        f  = -2*log(2)*s2*2^(-2*Bk) + lam*(2*alpha*Bk + beta);
        df =  4*log(2)^2*s2*2^(-2*Bk) + 2*lam*alpha;
        dB = f / (df + eps);
        Bk = max(B_floor, Bk - dB);
        if abs(dB) < 1e-11, break; end
    end
end

function verify_budgets(alloc_vecs, strategies, C_total, alpha)
    if nargin < 4, alpha = 1.37; end
    C_lut_u = alpha * sum(alloc_vecs{1}.^2);
    fprintf('  %-12s  %6s  %6s  %8s  %12s  %12s\n', ...
            'Strategy','min(B)','max(B)','sum(B)','sum(J(B))','C_lut_unif');
    for s = 1:numel(strategies)
        Bv = alloc_vecs{s};
        flags = '';
        if sum(Bv) ~= C_total, flags = [flags ' SUM_ERR']; end
        if any(Bv < 4),        flags = [flags ' FLOOR_VIOL']; end
        if any(Bv > 24),       flags = [flags ' CEIL_VIOL']; end
        if isempty(flags),     flags = ' OK'; end
        fprintf('  %-12s  %6d  %6d  %8d  %12.1f  %12.1f  %s\n', ...
                strategies{s}, min(Bv), max(Bv), sum(Bv), ...
                alpha*sum(Bv.^2), C_lut_u, flags);
    end
    fprintf('\n');
end