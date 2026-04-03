% e7_iso_budget - robustness across phonetically defined speech subsets.
%
% evaluates D_w for uniform and proposed across three subsets at N=835:
% whispered speech (bottom 15th percentile by rms), fricative-heavy
% (top 15th percentile by high-band energy above 4kHz), and unvoiced-burst
% (top 15th percentile by zero-crossing rate), plus the full corpus.
%
% also records per-clip acoustic features and per-band mse breakdowns
% to show which frequency regions drive any advantage or disadvantage.
%
% outputs:
%   E7_iso_summary.csv
%   E7_iso_acoustic_features.csv
%   E7_iso_band_breakdown.csv
%   E7_iso_clip_dist.csv
%
% run after E1 (requires E1_sigma2_eff.csv).
clear; clc; close all;

LIBRISPEECH_ROOTS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'
};

PERCENTILE = 15;        % raised from 5 -> N per subset ~835 instead of ~278
BUDGETS    = [40,60,80,100];

M=64; FS_TARGET=16000; B_UNIFORM=16; B_MAX=24; B_MIN=4;
ALPHA=1.37; BETA=0.00;
FRAME=2*M; HOP=M; WIN=hann(FRAME,'periodic');
C_TOTAL_REF=M*B_UNIFORM;
fc=(1:M)'*(FS_TARGET/2)/M;
S=exp(-0.5*(log(fc/2500)/0.55).^2); S=S/max(S);
high_bands=find(fc>4000);

% ── Load sigma2_eff from E1 (allocator input) ────────────────────────────
if ~isfile('E1_sigma2_eff.csv')
    error('E1_sigma2_eff.csv not found. Run E1 first.');
end
Te=readtable('E1_sigma2_eff.csv','VariableNamingRule','preserve');
sigma2_eff=Te.sigma2_eff;
fprintf('Loaded sigma2_eff from E1_sigma2_eff.csv.\n\n');

corpus=discover_librispeech(LIBRISPEECH_ROOTS); N_CLIPS=numel(corpus);

% ── Acoustic feature extraction ──────────────────────────────────────────
fprintf('Extracting acoustic features (%d clips)...\n', N_CLIPS);
rms_v=zeros(N_CLIPS,1); hbr_v=zeros(N_CLIPS,1);
zcr_v=zeros(N_CLIPS,1); valid=false(N_CLIPS,1);
for c=1:N_CLIPS
    try [x,fs]=audioread(corpus(c).path); catch; continue; end
    x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
    x=x/(max(abs(x))+eps);
    rms_v(c)=rms(x);
    Xk=wola_analyze(x,M,WIN,HOP);
    hbr_v(c)=sum(Xk(high_bands,:),'all')/(sum(Xk(:),'all')+eps);
    sg=sign(x); sg(sg==0)=1;
    zcr_v(c)=sum(abs(diff(sg)))/(2*length(x));
    valid(c)=true;
    if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
end
vidx=find(valid);
rv=rms_v(vidx); hv=hbr_v(vidx); zv=zcr_v(vidx);

% Subset masks using PERCENTILE
wm = rv <= prctile(rv, PERCENTILE);
fm = hv >= prctile(hv, 100-PERCENTILE);
um = zv >= prctile(zv, 100-PERCENTILE);

subsets  = {'Whispered','Fricative-Heavy','Unvoiced-Burst','Full-Corpus'};
masks    = {wm, fm, um, true(numel(vidx),1)};
strat_names = {'Uniform','Proposed'};
N_STRAT  = 2;   % Uniform and Proposed only; Linear omitted in results

fprintf('\nSubset sizes at PERCENTILE=%d:\n', PERCENTILE);
for ss=1:3
    fprintf('  %-22s  N=%d\n', subsets{ss}, sum(masks{ss}));
end
fprintf('  %-22s  N=%d\n\n', 'Full-Corpus', numel(vidx));

% Save acoustic features
fid=fopen('E7_iso_acoustic_features.csv','w');
fprintf(fid,'Clip,Speaker,Corpus,RMS,HighBandRatio,ZCR,IsWhispered,IsFricative,IsUnvoiced\n');
for i=1:numel(vidx); c=vidx(i);
    fprintf(fid,'%s,%s,%s,%.6f,%.6f,%.6f,%d,%d,%d\n', ...
            corpus(c).label,corpus(c).speaker,corpus(c).corpus, ...
            rv(i),hv(i),zv(i),wm(i),fm(i),um(i));
end
fclose(fid); fprintf('Saved: E7_iso_acoustic_features.csv\n\n');

% ── Pre-compute allocations ───────────────────────────────────────────────
fprintf('Computing allocations...\n');
N_BUD=numel(BUDGETS);
B_alloc=cell(N_BUD,N_STRAT);
for b=1:N_BUD
    pct=BUDGETS(b); C_total=round(pct/100*C_TOTAL_REF);
    B_alloc{b,1}=allocate_uniform(C_total,M);
    % Use Proposed (quadratic) for strategy 2
    B_alloc{b,2}=allocate_quadratic(sigma2_eff,ALPHA,BETA,C_total,M,B_MIN,B_MAX);
    verify_budgets(B_alloc(b,:),strat_names,C_total,ALPHA);
end

% ── Main evaluation loop ──────────────────────────────────────────────────
% Store: per-clip D_w [N_clips x N_STRAT] per budget x subset
% Store: per-band mean MSE [M x N_STRAT] per budget x subset
results    = struct();
clip_dist  = {};   % rows for E7_iso_clip_dist.csv
band_rows  = {};   % rows for E7_iso_band_breakdown.csv

for b=1:N_BUD
    pct=BUDGETS(b); C_total=round(pct/100*C_TOTAL_REF);
    fprintf('=== Budget %d%% (C=%d bits) ===\n', pct, C_total);

    for ss_idx=1:numel(subsets)
        ci_list = vidx(masks{ss_idx});
        N_sub   = numel(ci_list);
        dw_mat  = NaN(N_sub, N_STRAT);          % per-clip D_w
        band_sum= zeros(M, N_STRAT);             % sum of per-band MSE across clips
        band_cnt= 0;

        for ci=1:N_sub
            c=ci_list(ci);
            try [x,fs]=audioread(corpus(c).path);
            catch; continue; end
            x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
            x=x/(max(abs(x))+eps);
            Xk_ref=wola_analyze(x,M,WIN,HOP);

            for s=1:N_STRAT
                Bv=B_alloc{b,s};
                Xk_q=quant_sub(Xk_ref,Bv,B_MAX);
                band_err=zeros(M,1);
                for k=1:M
                    band_err(k)=mean((Xk_ref(k,:)-Xk_q(k,:)).^2);
                end
                dw_mat(ci,s)=sum(S.*band_err);
                band_sum(:,s)=band_sum(:,s)+band_err;
            end
            band_cnt=band_cnt+1;

            % Record clip-level data
            clip_dist{end+1}={pct, subsets{ss_idx}, corpus(c).label, ...
                              corpus(c).corpus, dw_mat(ci,1), dw_mat(ci,2)};
        end

        % Aggregate
        key=sprintf('s%d_b%d',ss_idx,pct);
        results.(key).subset  = subsets{ss_idx};
        results.(key).N       = N_sub;
        results.(key).med_dw  = median(dw_mat,'omitnan');
        results.(key).mean_dw = mean(dw_mat,'omitnan');
        results.(key).std_dw  = std(dw_mat,'omitnan');
        results.(key).ratio   = results.(key).med_dw(1) / ...
                                (results.(key).med_dw(2)+eps);
        results.(key).mean_band = band_sum / max(band_cnt,1);

        % Band breakdown rows
        for k=1:M
            band_rows{end+1}={pct, subsets{ss_idx}, k, fc(k), S(k), ...
                              results.(key).mean_band(k,1), ...
                              results.(key).mean_band(k,2)};
        end

        dB = 10*log10(results.(key).ratio);
        fprintf('  %-22s  N=%4d  median U/P=%.4f  (%+.2f dB)  mean U=%.3e P=%.3e\n', ...
                subsets{ss_idx}, N_sub, results.(key).ratio, dB, ...
                results.(key).mean_dw(1), results.(key).mean_dw(2));
    end
    fprintf('\n');
end

% ── Write E7_iso_summary.csv ──────────────────────────────────────────────
fid=fopen('E7_iso_summary.csv','w');
fprintf(fid,['BudgetPct,Subset,N,Uniform_Med_Dw,Proposed_Med_Dw,' ...
             'Uniform_Mean_Dw,Proposed_Mean_Dw,' ...
             'Uniform_Std_Dw,Proposed_Std_Dw,' ...
             'Uniform_Proposed_Ratio,dB\n']);
for b=1:N_BUD; pct=BUDGETS(b);
    for ss_idx=1:numel(subsets)
        r=results.(sprintf('s%d_b%d',ss_idx,pct));
        dB=10*log10(r.ratio);
        fprintf(fid,'%d,%s,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.4f,%.4f\n', ...
                pct,r.subset,r.N, ...
                r.med_dw(1),r.med_dw(2), ...
                r.mean_dw(1),r.mean_dw(2), ...
                r.std_dw(1),r.std_dw(2), ...
                r.ratio,dB);
    end
end
fclose(fid); fprintf('Saved: E7_iso_summary.csv\n');

% ── Write E7_iso_band_breakdown.csv ───────────────────────────────────────
fid=fopen('E7_iso_band_breakdown.csv','w');
fprintf(fid,'BudgetPct,Subset,BandIndex,FreqHz,SensS,Uniform_MeanMSE,Proposed_MeanMSE,Ratio_U_over_P,dB\n');
for i=1:numel(band_rows)
    r=band_rows{i};
    ratio_b=r{6}/(r{7}+eps);
    dB_b=10*log10(ratio_b);
    fprintf(fid,'%d,%s,%d,%.1f,%.6f,%.8e,%.8e,%.4f,%.4f\n', ...
            r{1},r{2},r{3},r{4},r{5},r{6},r{7},ratio_b,dB_b);
end
fclose(fid); fprintf('Saved: E7_iso_band_breakdown.csv\n');

% ── Write E7_iso_clip_dist.csv ────────────────────────────────────────────
fid=fopen('E7_iso_clip_dist.csv','w');
fprintf(fid,'BudgetPct,Subset,Clip,Corpus,Uniform_Dw,Proposed_Dw,Ratio_U_over_P\n');
for i=1:numel(clip_dist)
    r=clip_dist{i};
    rat=r{5}/(r{6}+eps);
    fprintf(fid,'%d,%s,%s,%s,%.8e,%.8e,%.4f\n', ...
            r{1},r{2},r{3},r{4},r{5},r{6},rat);
end
fclose(fid); fprintf('Saved: E7_iso_clip_dist.csv\n');

% ── Final console summary ─────────────────────────────────────────────────
fprintf('\n====== E7 FINAL SUMMARY (PERCENTILE=%d) ======\n\n', PERCENTILE);
fprintf('%-6s  %-22s  %5s  %8s  %8s\n', 'Budget','Subset','N','U/P','dB');
fprintf('%s\n', repmat('-',1,56));
for b=1:N_BUD; pct=BUDGETS(b);
    for ss_idx=1:numel(subsets)
        r=results.(sprintf('s%d_b%d',ss_idx,pct));
        dB=10*log10(r.ratio);
        fprintf('%-6d  %-22s  %5d  %8.4f  %+8.2f dB\n', ...
                pct, r.subset, r.N, r.ratio, dB);
    end
    fprintf('\n');
end

fprintf('Band breakdown saved to E7_iso_band_breakdown.csv\n');
fprintf('Inspect those bands with dB < 0 at 80%% budget to see\n');
fprintf('which frequency regions drive the disadvantage.\n');

% ── LOCAL FUNCTIONS ───────────────────────────────────────────────────────
function Xk=wola_analyze(x,M,win,hop)
    fl=2*M; nf=floor((length(x)-fl)/hop)+1; Xk=zeros(M,nf);
    for f=1:nf
        seg=x((f-1)*hop+1:(f-1)*hop+fl);
        sp=fft(seg.*win); Xk(:,f)=abs(sp(1:M)).^2;
    end
end

function Xk_q=quant_sub(Xk,Bv,Bmax)
    [M,n]=size(Xk); Xk_q=zeros(M,n);
    for k=1:M
        b=max(1,min(Bmax,round(Bv(k)))); sc=2^(b-1);
        Xk_q(k,:)=round(Xk(k,:)*sc)/sc;
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