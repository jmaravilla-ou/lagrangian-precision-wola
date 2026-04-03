% e9_iso_budget - haspi v2 and hasqi v2 evaluation.
%
% evaluates haspi v2 and hasqi v2 for uniform, linear, and proposed
% allocations at four budget levels on a stratified n=200 subsample.
% reports scores under both normal hearing and a mild-to-moderate
% sensorineural audiogram (25/30/35/45/55/60 dB HL at 250-6000 Hz).
%
% requires: haspi_v2.m, hasqi_v2.m (set haspi_dir below)
%           discover_librispeech.m on path
%
% outputs: E9_iso_haspi_<pct>pct.csv (x4), E9_iso_haspi_summary.csv
%
% run after E1 (requires E1_sigma2_eff.csv).

clear; clc; close all;

% haspi toolbox path - update this to point to your local copy
HASPI_DIR = '/path/to/haspi-toolbox';
if ~isfolder(HASPI_DIR)
    error('E9:path','HASPI folder not found:\n  %s', HASPI_DIR);
end
addpath(genpath(HASPI_DIR));
assert(exist('HASPI_v2','file')==2,'HASPI_v2.m not found on path.');
assert(exist('HASQI_v2','file')==2,'HASQI_v2.m not found on path.');

% parameters
M           = 64;
FS_TARGET   = 16000;
B_MAX       = 24;
B_MIN       = 4;
ALPHA       = 1.37;
BETA        = 0.00;
FRAME       = 2*M;     % 128 samples
HOP         = M;       % 64 samples
WIN         = hann(FRAME,'periodic');
C_TOTAL_REF = M * 16;  % 1024 bits
BUDGET_FRACS = [0.40, 0.60, 0.80, 1.00];

Level1   = 65;
HL_NH    = [0,  0,  0,  0,  0,  0];
HL_HI    = [25, 30, 35, 45, 55, 60];
EQ_FLAG  = 2;

SUBSAMPLE_N    = 200;
SUBSAMPLE_SEED = 42;

% psychoacoustic importance function - log-normal fit to ansi s3.5-1997 sii
% bin 1 = 125 hz, bin 64 = 8000 hz (no dc bin)
fc = (1:M)' * (FS_TARGET / 2) / M;
S  = exp(-0.5 * (log(fc / 2500) / 0.55).^2);
S  = S / max(S);

% corpus
CORPUS_DIRS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'};
corpus = discover_librispeech(CORPUS_DIRS);

% stratified subsample
rng(SUBSAMPLE_SEED);
idx        = randperm(numel(corpus), min(SUBSAMPLE_N, numel(corpus)));
corpus_sub = corpus(idx);
N_SUB      = numel(corpus_sub);
fprintf('Subsampled %d clips for HASPI/HASQI evaluation.\n\n', N_SUB);

% long-term subband power estimate
fprintf('Computing LTSS across %d clips...\n', numel(corpus));
sigma2_sum = zeros(M,1); n_valid = 0;
for c = 1:numel(corpus)
    try [x,fs] = audioread(corpus(c).path); catch; continue; end
    x = mean(x,2);
    if fs ~= FS_TARGET, x = resample(x,FS_TARGET,fs); end
    x = x / (max(abs(x))+eps);
    Xk = wola_analyze(x,M,WIN,HOP);
    sigma2_sum = sigma2_sum + mean(Xk,2);
    n_valid = n_valid + 1;
    if mod(c,500)==0, fprintf('  LTSS: %d/%d\n',c,numel(corpus)); end
end
raw        = sigma2_sum/n_valid + eps;
raw        = raw / max(raw);
sigma2_eff = raw .* S;
sigma2_eff = sigma2_eff / max(sigma2_eff);
fprintf('LTSS done.\n\n');

% pre-compute allocations
strategies = {'Uniform','Linear','Proposed'};
N_STRAT    = 3;
n_bud      = numel(BUDGET_FRACS);
B_alloc    = cell(n_bud, N_STRAT);

for b = 1:n_bud
    C_total        = round(C_TOTAL_REF * BUDGET_FRACS(b));
    B_alloc{b,1}   = allocate_uniform(C_total, M);
    B_alloc{b,2}   = allocate_linear(sigma2_eff, C_total, M, B_MIN, B_MAX);
    B_alloc{b,3}   = allocate_quadratic(sigma2_eff, ALPHA, BETA, C_total, M, B_MIN, B_MAX);
    fprintf('Budget %3.0f%%:', BUDGET_FRACS(b)*100);
    for s = 1:N_STRAT
        fprintf('  %s sum=%d', strategies{s}, sum(B_alloc{b,s}));
    end
    fprintf('  (target %d)\n', C_total);
    verify_budgets(B_alloc(b,:), strategies, C_total, ALPHA);
end

% main evaluation loop
% results: [N_SUB x 8 x n_bud]
% cols: HASPI_NH_U, HASPI_NH_L, HASPI_NH_P,
%       HASPI_HI_U, HASPI_HI_L, HASPI_HI_P,
%       HASQI_NH_U, HASQI_NH_L, HASQI_NH_P,
%       HASQI_HI_U, HASQI_HI_L, HASQI_HI_P
N_METRICS = 12;
results   = NaN(N_SUB, N_METRICS, n_bud);
t_start   = tic;

for b = 1:n_bud
    pct = BUDGET_FRACS(b) * 100;
    fprintf('=== Budget %.0f%% ===\n', pct);

    for c = 1:N_SUB
        try
            [x_raw,fs] = audioread(corpus_sub(c).path);
        catch me
            warning('Clip %d skipped: %s', c, me.message); continue;
        end
        x_raw = mean(x_raw,2);
        if fs ~= FS_TARGET, x_raw = resample(x_raw,FS_TARGET,fs); end
        x_ref = x_raw / (rms(x_raw)+eps);

        [Xk_mag, ph] = wola_analyze_complex(x_ref, M, WIN, HOP);
        sig_len       = length(x_ref);

        for s = 1:N_STRAT
            Bv   = B_alloc{b,s};
            Xk_q = quantize_subbands(Xk_mag, Bv, B_MAX);
            y    = wola_synthesize(Xk_q, ph, M, WIN, HOP, sig_len);
            y    = max(-1, min(1, y));
            y    = y / (rms(y)+eps);

            col_off = s - 1;   % 0-based offset within each metric group

            try [hi_nh,~] = HASPI_v2(x_ref,FS_TARGET,y,FS_TARGET,HL_NH,Level1);
            catch; hi_nh = NaN; end
            try [hi_hi,~] = HASPI_v2(x_ref,FS_TARGET,y,FS_TARGET,HL_HI,Level1);
            catch; hi_hi = NaN; end
            try [qi_nh,~,~,~] = HASQI_v2(x_ref,FS_TARGET,y,FS_TARGET,HL_NH,EQ_FLAG,Level1);
            catch; qi_nh = NaN; end
            try [qi_hi,~,~,~] = HASQI_v2(x_ref,FS_TARGET,y,FS_TARGET,HL_HI,EQ_FLAG,Level1);
            catch; qi_hi = NaN; end

            results(c, 1+col_off, b) = hi_nh;
            results(c, 4+col_off, b) = hi_hi;
            results(c, 7+col_off, b) = qi_nh;
            results(c,10+col_off, b) = qi_hi;
        end

        if mod(c,50)==0 || c==N_SUB
            rm  = nanmean(results(1:c,:,b), 1);
            ela = toc(t_start);
            eta = ela/((b-1)*N_SUB+c) * (n_bud*N_SUB-(b-1)*N_SUB-c);
            fprintf('  %3d/%d  HASPI NH U/L/P=%.3f/%.3f/%.3f  HI=%.3f/%.3f/%.3f  [ETA %.0fs]\n', ...
                c, N_SUB, rm(1),rm(2),rm(3), rm(4),rm(5),rm(6), eta);
        end
    end

    % Per-budget CSV
    pct_str = sprintf('%dpct', round(pct));
    fname   = sprintf('E9_iso_haspi_%s.csv', pct_str);
    fid     = fopen(fname,'w');
    fprintf(fid,['Clip,Speaker,Corpus,' ...
        'HASPI_NH_Uniform,HASPI_NH_Linear,HASPI_NH_Proposed,' ...
        'HASPI_HI_Uniform,HASPI_HI_Linear,HASPI_HI_Proposed,' ...
        'HASQI_NH_Uniform,HASQI_NH_Linear,HASQI_NH_Proposed,' ...
        'HASQI_HI_Uniform,HASQI_HI_Linear,HASQI_HI_Proposed\n']);
    for c = 1:N_SUB
        [~,cname,~] = fileparts(corpus_sub(c).path);
        fprintf(fid,'%s,%s,%s', cname, corpus_sub(c).speaker, corpus_sub(c).corpus);
        for j = 1:N_METRICS, fprintf(fid,',%.6f',results(c,j,b)); end
        fprintf(fid,'\n');
    end
    fclose(fid);
    fprintf('  Saved %s\n\n', fname);
end

% summary
fid_s = fopen('E9_iso_haspi_summary.csv','w');
fprintf(fid_s,'BudgetPct,Audiogram,Strategy,HASPI_Mean,HASPI_Std,HASQI_Mean,HASQI_Std\n');
labels = {'NH','HI'};
pi_cols = {[1 2 3],[4 5 6]};
qi_cols = {[7 8 9],[10 11 12]};

fprintf('\n======= HASPI/HASQI Summary =======\n');
for b = 1:n_bud
    pct = BUDGET_FRACS(b)*100;
    R   = results(:,:,b);
    for a = 1:2
        for s = 1:N_STRAT
            pi_v = R(:, pi_cols{a}(s));
            qi_v = R(:, qi_cols{a}(s));
            fprintf(fid_s,'%.0f,%s,%s,%.6f,%.6f,%.6f,%.6f\n', ...
                pct, labels{a}, strategies{s}, ...
                nanmean(pi_v), nanstd(pi_v), ...
                nanmean(qi_v), nanstd(qi_v));
        end
    end
    mu = nanmean(R,1);
    fprintf('Budget %3.0f%%  |  HASPI HI:  U=%.4f  L=%.4f  P=%.4f\n', ...
        pct, mu(4), mu(5), mu(6));
    fprintf('             |  HASQI HI:  U=%.4f  L=%.4f  P=%.4f\n', ...
        mu(10), mu(11), mu(12));
end
fclose(fid_s);
fprintf('Saved: E9_iso_haspi_summary.csv\n');
fprintf('Done. Elapsed: %.1f s\n', toc(t_start));

% local functions
function Xk = wola_analyze(x,M,win,hop)
    fl=2*M; nf=floor((length(x)-fl)/hop)+1; Xk=zeros(M,nf);
    for f=1:nf
        seg=x((f-1)*hop+1:(f-1)*hop+fl);
        sp=fft(seg.*win); Xk(:,f)=abs(sp(1:M)).^2;
    end
end

function [Xk_mag,phase] = wola_analyze_complex(x,M,win,hop)
    fl=2*M; nf=floor((length(x)-fl)/hop)+1;
    Xk_mag=zeros(M,nf); phase=zeros(M,nf);
    for f=1:nf
        seg=x((f-1)*hop+1:(f-1)*hop+fl);
        sp=fft(seg.*win); sh=sp(1:M);
        Xk_mag(:,f)=abs(sh); phase(:,f)=angle(sh);
    end
end

function x_out = wola_synthesize(Xk_mag,phase,M,win,hop,sig_len)
    fl=2*M; nf=size(Xk_mag,2);
    x_out=zeros(sig_len+fl,1); na=zeros(sig_len+fl,1);
    for f=1:nf
        sh=Xk_mag(:,f).*exp(1j*phase(:,f));
        sf=zeros(fl,1);
        sf(1)=real(sh(1)); sf(2:M)=sh(2:M); sf(M+1)=0;
        sf(M+2:end)=conj(flipud(sh(2:M)));
        fr=real(ifft(sf)).*win;
        idx=(f-1)*hop+(1:fl); idx=idx(idx<=length(x_out));
        x_out(idx)=x_out(idx)+fr(1:length(idx));
        na(idx)=na(idx)+win(1:length(idx)).^2;
    end
    na(na<1e-8)=1; x_out=x_out(1:sig_len)./na(1:sig_len);
end

function Xk_q = quantize_subbands(Xk_ref,Bv,B_max)
    [M,nf]=size(Xk_ref); Xk_q=zeros(M,nf);
    for k=1:M
        bk=max(1,min(B_max,round(Bv(k)))); sc=2^(bk-1);
        Xk_q(k,:)=round(Xk_ref(k,:)*sc)/sc;
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