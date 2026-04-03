% e11_bitdepth_sweep - per-band bit-depth sensitivity and B_min derivation.
%
% for each of the 64 subbands independently, sweeps bit depth from B=24
% down to B=1. at each step only that band is quantized; all others stay
% at full 24-bit precision. a companion global sweep quantizes all bands.
%
% the original input x is used as the reference throughout, not the wola
% reconstruction at B=24, so snr and stoi reflect absolute quality.
% phase is stored from analysis and used in synthesis.
%
% metrics: snr (db), stoi, per-band subband mse, per-band snr vs theory,
% and pesq (if matlab audio toolbox is available).
%
% outputs:
%   E11_perband_snr.csv
%   E11_perband_stoi.csv
%   E11_perband_reconerr.csv
%   E11_perband_freqdev.csv
%   E11_perband_pesq.csv  (if pesq available)
%   E11_global_sweep.csv
%   E11_bmin_summary.csv
%   E11_summary_report.txt
clear; clc; close all;

LIBRISPEECH_ROOTS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'
};

% ── Parameters ────────────────────────────────────────────────────────────
M              = 64;
FS_TARGET      = 16000;
B_MAX          = 24;
FRAME          = 2*M;
HOP            = M;
WIN            = hann(FRAME, 'periodic');
B_SWEEP        = 1:B_MAX;
N_SWEEP        = numel(B_SWEEP);
N_CLIPS_TARGET = 100;

% Quality thresholds for B_min derivation
SNR_THRESH     = 20.0;   % dB  — noise clearly audible below this
STOI_THRESH    = 0.95;   % —— intelligibility begins to drop below this
PESQ_THRESH    = 3.0;    % MOS — "fair" boundary on 1–4.5 scale

fc = (1:M)' * (FS_TARGET/2) / M;
S  = exp(-0.5*(log(fc/2500)/0.55).^2);  S = S/max(S);

% ── Toolbox availability ───────────────────────────────────────────────────
USE_STOI = exist('stoi', 'file') == 2;
USE_PESQ = exist('pesq', 'file') == 2;
fprintf('STOI available: %d\n', USE_STOI);
fprintf('PESQ available: %d\n\n', USE_PESQ);
if ~USE_STOI
    fprintf('WARNING: stoi() not found. STOI column will be zeros.\n');
    fprintf('         Install the Perceptual Audio Evaluation toolbox or\n');
    fprintf('         add stoi.m to your MATLAB path.\n\n');
end

% ── Corpus: stratified subsample ──────────────────────────────────────────
corpus_full = discover_librispeech(LIBRISPEECH_ROOTS);
N_FULL = numel(corpus_full);
step   = max(1, floor(N_FULL / N_CLIPS_TARGET));
idx    = 1:step:N_FULL;
idx    = idx(1:min(N_CLIPS_TARGET, numel(idx)));
N_CLIPS = numel(idx);
corpus  = corpus_full(idx);
fprintf('Corpus: %d clips sampled from %d (step=%d).\n\n', N_CLIPS, N_FULL, step);

% ── Pre-load clips and run WOLA analysis once per clip ────────────────────
% Store FULL COMPLEX subband matrix so phase is available for synthesis.
% Also store x for use as the SNR/STOI/PESQ reference.
fprintf('Pre-loading and analysing %d clips...\n', N_CLIPS);
clips    = cell(N_CLIPS, 1);   % original waveform x
Xk_full  = cell(N_CLIPS, 1);   % [M x N_frames] COMPLEX subband matrix

n_ok = 0;
for c = 1:N_CLIPS
    try
        [x, fs] = audioread(corpus(c).path);
        x = mean(x, 2);
        if fs ~= FS_TARGET
            x = resample(x, FS_TARGET, fs);
        end
        x = x / (max(abs(x)) + eps);
        clips{c}   = x;
        Xk_full{c} = wola_analyze(x, M, WIN, HOP);   % complex [M x N_frames]
        n_ok = n_ok + 1;
    catch e
        clips{c}   = [];
        Xk_full{c} = [];
        fprintf('  WARNING: failed to load clip %d: %s\n', c, e.message);
    end
    if mod(c, 25) == 0
        fprintf('  Loaded %d/%d\n', c, N_CLIPS);
    end
end
fprintf('Pre-load complete: %d/%d clips valid.\n\n', n_ok, N_CLIPS);

% Quick sanity check: reconstruct clip 1 at full precision and measure SNR
% This confirms the WOLA round-trip fidelity before the sweep begins.
fprintf('--- WOLA round-trip sanity check (clip 1, B=24) ---\n');
c_test = find(~cellfun(@isempty, clips), 1);
if ~isempty(c_test)
    x_t  = clips{c_test};
    Xk_t = Xk_full{c_test};
    y_t  = wola_synthesize(Xk_t, M, WIN, HOP, length(x_t));
    rt_snr = 10*log10(mean(x_t.^2) / (mean((x_t - y_t).^2) + eps));
    fprintf('  Round-trip SNR vs original x: %.2f dB\n', rt_snr);
    fprintf('  (Expected: ~30-60 dB for a correct WOLA implementation)\n');
    fprintf('  (If this is negative, WOLA synthesis has a bug — stop here)\n\n');
    if rt_snr < 0
        error(['WOLA round-trip SNR is %.2f dB — synthesis is broken. ' ...
               'Check wola_analyze and wola_synthesize implementations.'], rt_snr);
    end
else
    error('No valid clips loaded. Check LIBRISPEECH_ROOTS paths.');
end

% ── Accumulators ──────────────────────────────────────────────────────────
% [M x N_SWEEP] matrices, one entry per (band, bit-depth) condition
pb_snr_sum     = zeros(M, N_SWEEP);
pb_stoi_sum    = zeros(M, N_SWEEP);
pb_err_sum     = zeros(M, N_SWEEP);
pb_fdev_sum    = zeros(M, N_SWEEP);
pb_pesq_sum    = zeros(M, N_SWEEP);
pb_cnt         = zeros(M, N_SWEEP);   % clips contributing to SNR/STOI/err
pb_pesq_cnt    = zeros(M, N_SWEEP);   % clips contributing to PESQ
pb_stoi_cnt    = zeros(M, N_SWEEP);   % clips contributing to STOI

% ── PART 1: PER-BAND SWEEP ────────────────────────────────────────────────
fprintf('=== PART 1: Per-Band Bit-Depth Sweep ===\n');
fprintf('%d bands x %d bit depths x %d clips = %d conditions\n\n', ...
        M, N_SWEEP, N_CLIPS, M*N_SWEEP);

total = M * N_SWEEP;
done  = 0;
t0    = tic;

for k = 1:M
    for bi = 1:N_SWEEP
        B  = B_SWEEP(bi);
        sc = 2^(B - 1);

        for c = 1:N_CLIPS
            if isempty(clips{c}), continue; end

            x  = clips{c};
            Xk = Xk_full{c};        % full complex [M x N_frames]
            Lx = length(x);

            % ── Quantize only band k ──────────────────────────────────
            % Copy Xk, replace row k with quantized magnitude + original phase
            Xk_q        = Xk;
            mag_k       = abs(Xk(k,:));
            phase_k     = angle(Xk(k,:));
            mag_k_q     = round(mag_k * sc) / sc;
            Xk_q(k,:)   = mag_k_q .* exp(1j * phase_k);

            % ── Reconstruct ──────────────────────────────────────────
            y_q = wola_synthesize(Xk_q, M, WIN, HOP, Lx);

            % ── Metric 1: Output SNR vs original x ───────────────────
            sp  = mean(x.^2) + eps;
            ep  = mean((x - y_q).^2) + eps;
            pb_snr_sum(k,bi) = pb_snr_sum(k,bi) + 10*log10(sp/ep);

            % ── Metric 2: STOI vs original x ─────────────────────────
            if USE_STOI && Lx >= FS_TARGET * 0.3
                try
                    sv = stoi(x, y_q, FS_TARGET);
                    if isfinite(sv)
                        pb_stoi_sum(k,bi)  = pb_stoi_sum(k,bi) + sv;
                        pb_stoi_cnt(k,bi)  = pb_stoi_cnt(k,bi) + 1;
                    end
                catch
                    % skip — do not inject a fake value
                end
            end

            % ── Metric 3: Per-band subband MSE ───────────────────────
            pb_err_sum(k,bi) = pb_err_sum(k,bi) + ...
                               mean((mag_k_q - mag_k).^2);

            % ── Metric 4: Per-band SNR vs theoretical ────────────────
            sig_k = mean(mag_k.^2) + eps;
            err_k = mean((mag_k_q - mag_k).^2) + eps;
            snr_actual_k  = 10*log10(sig_k / err_k);
            snr_theory_k  = 6.02*B + 1.76;
            pb_fdev_sum(k,bi) = pb_fdev_sum(k,bi) + ...
                                abs(snr_theory_k - snr_actual_k);

            % ── Metric 5: PESQ vs original x ─────────────────────────
            if USE_PESQ && Lx >= FS_TARGET * 0.5
                try
                    mv = pesq(FS_TARGET, x, y_q, 'wb');
                    if isfinite(mv)
                        pb_pesq_sum(k,bi)  = pb_pesq_sum(k,bi) + mv;
                        pb_pesq_cnt(k,bi)  = pb_pesq_cnt(k,bi) + 1;
                    end
                catch
                    % skip — do not inject a fake value
                end
            end

            pb_cnt(k,bi) = pb_cnt(k,bi) + 1;
        end  % clips

        done = done + 1;
        if mod(done, 64) == 0
            el  = toc(t0);
            eta = (total - done) / (done/el);
            fprintf('  %d/%d conditions  elapsed=%.0fs  ETA=%.0fs\n', ...
                    done, total, el, eta);
        end
    end  % bit depths
end  % bands

fprintf('Part 1 done (%.1f s).\n\n', toc(t0));

% ── Average accumulators into mean matrices ───────────────────────────────
pb_snr  = pb_snr_sum  ./ max(pb_cnt, 1);
pb_err  = pb_err_sum  ./ max(pb_cnt, 1);
pb_fdev = pb_fdev_sum ./ max(pb_cnt, 1);

pb_stoi = zeros(M, N_SWEEP);
for k = 1:M
    for bi = 1:N_SWEEP
        if pb_stoi_cnt(k,bi) > 0
            pb_stoi(k,bi) = pb_stoi_sum(k,bi) / pb_stoi_cnt(k,bi);
        else
            pb_stoi(k,bi) = NaN;
        end
    end
end

pb_pesq = zeros(M, N_SWEEP);
for k = 1:M
    for bi = 1:N_SWEEP
        if pb_pesq_cnt(k,bi) > 0
            pb_pesq(k,bi) = pb_pesq_sum(k,bi) / pb_pesq_cnt(k,bi);
        else
            pb_pesq(k,bi) = NaN;
        end
    end
end

% ── PART 2: GLOBAL UNIFORM SWEEP ──────────────────────────────────────────
fprintf('=== PART 2: Global Uniform Sweep ===\n');

glob_snr      = zeros(N_SWEEP, 1);
glob_stoi     = zeros(N_SWEEP, 1);
glob_pesq     = zeros(N_SWEEP, 1);
glob_cnt      = zeros(N_SWEEP, 1);
glob_stoi_cnt = zeros(N_SWEEP, 1);
glob_pesq_cnt = zeros(N_SWEEP, 1);

for bi = 1:N_SWEEP
    B  = B_SWEEP(bi);
    sc = 2^(B - 1);

    snr_sum = 0;  stoi_sum = 0;  pesq_sum = 0;
    cnt = 0;  stoi_cnt = 0;  pesq_cnt = 0;

    for c = 1:N_CLIPS
        if isempty(clips{c}), continue; end
        x  = clips{c};
        Xk = Xk_full{c};
        Lx = length(x);

        % Quantize ALL bands to B bits
        mag_q = round(abs(Xk) * sc) / sc;
        Xk_q  = mag_q .* exp(1j * angle(Xk));
        y_q   = wola_synthesize(Xk_q, M, WIN, HOP, Lx);

        sp = mean(x.^2) + eps;
        ep = mean((x - y_q).^2) + eps;
        snr_sum = snr_sum + 10*log10(sp/ep);
        cnt     = cnt + 1;

        if USE_STOI && Lx >= FS_TARGET*0.3
            try
                sv = stoi(x, y_q, FS_TARGET);
                if isfinite(sv), stoi_sum = stoi_sum+sv; stoi_cnt = stoi_cnt+1; end
            catch; end
        end

        if USE_PESQ && Lx >= FS_TARGET*0.5
            try
                mv = pesq(FS_TARGET, x, y_q, 'wb');
                if isfinite(mv), pesq_sum = pesq_sum+mv; pesq_cnt = pesq_cnt+1; end
            catch; end
        end
    end

    if cnt > 0
        glob_snr(bi)  = snr_sum / cnt;
        glob_cnt(bi)  = cnt;
    end
    if stoi_cnt > 0
        glob_stoi(bi)     = stoi_sum / stoi_cnt;
        glob_stoi_cnt(bi) = stoi_cnt;
    else
        glob_stoi(bi) = NaN;
    end
    if pesq_cnt > 0
        glob_pesq(bi)     = pesq_sum / pesq_cnt;
        glob_pesq_cnt(bi) = pesq_cnt;
    else
        glob_pesq(bi) = NaN;
    end

    fprintf('  B=%2d: SNR=%7.2f dB  STOI=%.4f  PESQ=%.3f\n', ...
            B, glob_snr(bi), glob_stoi(bi), glob_pesq(bi));
end
fprintf('Part 2 done.\n\n');

% ── PART 3: B_min DERIVATION ──────────────────────────────────────────────
fprintf('=== PART 3: B_min Derivation ===\n\n');

bmin_snr    = zeros(M,1);
bmin_stoi   = zeros(M,1);
bmin_pesq   = zeros(M,1);
bmin_global = zeros(M,1);

for k = 1:M
    bmin_snr(k)  = find_bmin(pb_snr(k,:),  B_SWEEP, SNR_THRESH,  'above');
    bmin_stoi(k) = find_bmin(pb_stoi(k,:), B_SWEEP, STOI_THRESH, 'above');
    bmin_pesq(k) = find_bmin(pb_pesq(k,:), B_SWEEP, PESQ_THRESH, 'above');
    % Global = most conservative of available metrics
    candidates = [bmin_snr(k)];
    if USE_STOI, candidates(end+1) = bmin_stoi(k); end %#ok<AGROW>
    if USE_PESQ, candidates(end+1) = bmin_pesq(k); end %#ok<AGROW>
    bmin_global(k) = max(candidates);
end

% Print summary table
fprintf('%-6s %-10s %-8s %-10s %-10s %-10s %-12s\n', ...
        'Band','Freq(Hz)','S(k)','Bmin_SNR','Bmin_STOI','Bmin_PESQ','Bmin_Global');
fprintf('%s\n', repmat('-',1,68));
for k = 1:M
    fprintf('k=%-3d  %6.0f Hz  %.3f  %8d  %9d  %9d  %10d\n', ...
            k, fc(k), S(k), bmin_snr(k), bmin_stoi(k), bmin_pesq(k), bmin_global(k));
end
sys_bmin = max(bmin_global);
fprintf('\nSystem-wide B_min: %d bits\n\n', sys_bmin);

% ── WRITE OUTPUT FILES ────────────────────────────────────────────────────
fprintf('Writing output files...\n');

write_matrix_csv('E11_perband_snr.csv',      pb_snr,  M, B_SWEEP, 'SNR_dB',    fc);
write_matrix_csv('E11_perband_stoi.csv',     pb_stoi, M, B_SWEEP, 'STOI',      fc);
write_matrix_csv('E11_perband_reconerr.csv', pb_err,  M, B_SWEEP, 'ReconErr',  fc);
write_matrix_csv('E11_perband_freqdev.csv',  pb_fdev, M, B_SWEEP, 'FreqDev_dB',fc);
if USE_PESQ
    write_matrix_csv('E11_perband_pesq.csv', pb_pesq, M, B_SWEEP, 'PESQ_MOS',  fc);
end

fid = fopen('E11_global_sweep.csv','w');
fprintf(fid,'BitDepth,SNR_dB,STOI,PESQ_MOS,N_clips\n');
for bi = 1:N_SWEEP
    fprintf(fid,'%d,%.6f,%.6f,%.6f,%d\n', ...
            B_SWEEP(bi), glob_snr(bi), glob_stoi(bi), glob_pesq(bi), glob_cnt(bi));
end
fclose(fid);
fprintf('  Saved: E11_global_sweep.csv\n');

fid = fopen('E11_bmin_summary.csv','w');
fprintf(fid,'BandIndex,FreqHz,SensS,Bmin_SNR,Bmin_STOI,Bmin_PESQ,Bmin_Global\n');
for k = 1:M
    fprintf(fid,'%d,%.1f,%.6f,%d,%d,%d,%d\n', ...
            k, fc(k), S(k), bmin_snr(k), bmin_stoi(k), bmin_pesq(k), bmin_global(k));
end
fclose(fid);
fprintf('  Saved: E11_bmin_summary.csv\n');

fid = fopen('E11_summary_report.txt','w');
fprintf(fid,'E11 BIT-DEPTH SWEEP — SUMMARY REPORT\n');
fprintf(fid,'%s\n\n', repmat('=',1,50));
fprintf(fid,'WOLA round-trip SNR (B=24 vs original x): %.2f dB\n\n', rt_snr);
fprintf(fid,'GLOBAL UNIFORM SWEEP:\n');
fprintf(fid,'%-8s %-12s %-10s %-10s\n','B(bits)','SNR(dB)','STOI','PESQ');
for bi = 1:N_SWEEP
    fprintf(fid,'%-8d %-12.2f %-10.4f %-10.3f\n', ...
            B_SWEEP(bi), glob_snr(bi), glob_stoi(bi), glob_pesq(bi));
end
fprintf(fid,'\nPER-BAND B_min SUMMARY:\n');
fprintf(fid,'%-6s %-10s %-8s %-10s %-10s %-10s %-12s\n', ...
        'Band','Freq(Hz)','S(k)','Bmin_SNR','Bmin_STOI','Bmin_PESQ','Bmin_Global');
for k = 1:M
    fprintf(fid,'%-6d %-10.0f %-8.4f %-10d %-10d %-10d %-12d\n', ...
            k, fc(k), S(k), bmin_snr(k), bmin_stoi(k), bmin_pesq(k), bmin_global(k));
end
fprintf(fid,'\nSYSTEM-WIDE B_min: %d bits\n', sys_bmin);
fclose(fid);
fprintf('  Saved: E11_summary_report.txt\n\n');

fprintf('=== E11 COMPLETE ===\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function Xk = wola_analyze(x, M, win, hop)
    % Returns FULL COMPLEX [M x N_frames] subband matrix.
    % Phase is preserved — do NOT discard it.
    fl = 2*M;
    nf = floor((length(x) - fl) / hop) + 1;
    Xk = zeros(M, nf);   % complex
    for f = 1:nf
        i0  = (f-1)*hop + 1;
        seg = x(i0 : i0+fl-1);
        sp  = fft(seg .* win);
        Xk(:,f) = sp(1:M);   % complex — magnitude AND phase
    end
end

function y = wola_synthesize(Xk, M, win, hop, slen)
    % Reconstructs time-domain signal from complex subband matrix Xk.
    % Xk: [M x N_frames] complex (full spectrum stored)
    fl   = 2*M;
    nf   = size(Xk, 2);
    buf  = zeros(slen + fl, 1);
    norm = zeros(slen + fl, 1);
    for f = 1:nf
        % Rebuild full symmetric spectrum for real-valued IFFT
        sh  = Xk(:, f);              % M complex bins
        sf  = zeros(fl, 1);
        sf(1)         = real(sh(1));       % DC — real
        sf(2:M)       = sh(2:M);           % positive frequencies
        sf(M+1)       = 0;                 % Nyquist — zero
        sf(M+2:end)   = conj(flipud(sh(2:M)));  % conjugate-symmetric
        fr  = real(ifft(sf)) .* win;
        i0  = (f-1)*hop + 1;
        i1  = min(i0+fl-1, length(buf));
        n   = i1 - i0 + 1;
        buf(i0:i1)  = buf(i0:i1)  + fr(1:n);
        norm(i0:i1) = norm(i0:i1) + win(1:n).^2;
    end
    norm(norm < 1e-8) = 1;
    y = buf(1:slen) ./ norm(1:slen);
end

function bmin = find_bmin(metric_vec, B_sweep, threshold, direction)
    % Scan B_sweep(1)..B_sweep(end) and return the FIRST B that satisfies
    % the threshold. For 'above': first B where metric >= threshold.
    % Returns B_sweep(end) if no B satisfies the condition.
    bmin = B_sweep(end);
    for i = 1:numel(B_sweep)
        v = metric_vec(i);
        if isnan(v), continue; end
        if strcmp(direction,'above') && v >= threshold
            bmin = B_sweep(i); return;
        elseif strcmp(direction,'below') && v <= threshold
            bmin = B_sweep(i); return;
        end
    end
end

function write_matrix_csv(fname, mat, M, B_sweep, metric_name, fc)
    fid = fopen(fname, 'w');
    fprintf(fid, 'BandIndex,FreqHz');
    for bi = 1:numel(B_sweep)
        fprintf(fid, ',%s_B%d', metric_name, B_sweep(bi));
    end
    fprintf(fid, '\n');
    for k = 1:M
        fprintf(fid, '%d,%.1f', k, fc(k));
        for bi = 1:numel(B_sweep)
            v = mat(k,bi);
            if isnan(v)
                fprintf(fid, ',NaN');
            else
                fprintf(fid, ',%.8f', v);
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  Saved: %s\n', fname);
end