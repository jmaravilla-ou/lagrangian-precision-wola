% e2_iso_budget - architectural parity isolation.
%
% confirms the D_w advantage is purely algorithmic by measuring
% mid-band and high-band envelope distortion per strategy at each budget.
%
% outputs: E2_iso_summary.csv, E2_iso_alloc_<pct>pct.csv (x4)
clear; clc; close all;

LIBRISPEECH_ROOTS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'
};

M=64; FS_TARGET=16000; B_UNIFORM=16; B_MAX=24; B_MIN=4;
ALPHA=1.37; BETA=0.00;
FRAME=2*M; HOP=M; WIN=hann(FRAME,'periodic');
C_TOTAL_REF=M*B_UNIFORM;
BUDGET_FRACS=[0.40,0.60,0.80,1.00];
fc=(1:M)'*(FS_TARGET/2)/M;
S=exp(-0.5*(log(fc/2500)/0.55).^2); S=S/max(S);
mid_bands=find(fc>=1000&fc<=4000); high_bands=find(fc>4000);

strategies={'Uniform','Linear','Proposed'};
N_STRAT=3; N_BUD=numel(BUDGET_FRACS);

corpus=discover_librispeech(LIBRISPEECH_ROOTS); N_CLIPS=numel(corpus);

fprintf('Computing LTSS...\n');
ss=zeros(M,1); nv=0;
for c=1:N_CLIPS
    try [x,fs]=audioread(corpus(c).path); catch; continue; end
    x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
    x=x/(max(abs(x))+eps); Xk=wola_analyze(x,M,WIN,HOP);
    ss=ss+mean(Xk,2); nv=nv+1;
    if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
end
raw=ss/nv+eps; raw=raw/max(raw);
sigma2_eff=raw.*S; sigma2_eff=sigma2_eff/max(sigma2_eff);
fprintf('LTSS done.\n\n');

summary_data=cell(N_BUD,N_STRAT);

for b=1:N_BUD
    frac=BUDGET_FRACS(b); pct=round(frac*100); C_total=round(frac*C_TOTAL_REF);
    fprintf('[%d%%] C=%d bits\n',pct,C_total);

    B_u=allocate_uniform(C_total,M);
    B_l=allocate_linear(sigma2_eff,C_total,M,B_MIN,B_MAX);
    B_q=allocate_quadratic(sigma2_eff,ALPHA,BETA,C_total,M,B_MIN,B_MAX);
    alloc_vecs={B_u,B_l,B_q};

    % Verification
    verify_budgets(alloc_vecs, strategies, C_total, ALPHA);

    fname=sprintf('E2_iso_alloc_%dpct.csv',pct);
    fid=fopen(fname,'w');
    fprintf(fid,'BandIndex,FreqHz,Uniform,Linear,Proposed\n');
    for k=1:M
        fprintf(fid,'%d,%.1f,%.0f,%.0f,%.0f\n',k,fc(k),B_u(k),B_l(k),B_q(k));
    end
    fclose(fid);

    corr=zeros(N_CLIPS,N_STRAT,2);
    for c=1:N_CLIPS
        try [x,fs]=audioread(corpus(c).path); catch; corr(c,:,:)=NaN; continue; end
        x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
        x=x/(max(abs(x))+eps); Xk_ref=wola_analyze(x,M,WIN,HOP);
        for s=1:N_STRAT
            Bv=alloc_vecs{s}; Xk_q=quant_sub(Xk_ref,Bv,B_MAX);
            corr(c,s,1)=env_dist(rms(Xk_ref(mid_bands,:),2), ...
                                  rms(Xk_q(mid_bands,:),2));
            corr(c,s,2)=env_dist(rms(Xk_ref(high_bands,:),2), ...
                                  rms(Xk_q(high_bands,:),2));
        end
        if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
    end

    for s=1:N_STRAT
        Bv=alloc_vecs{s}; mc=corr(:,s,1); hc=corr(:,s,2);
        summary_data{b,s}=struct('budget_pct',pct,'strategy',strategies{s}, ...
            'avg_bits',mean(Bv),'high_bits',mean(Bv(high_bands)), ...
            'mid_mean',mean(mc,'omitnan'),'mid_std',std(mc,'omitnan'), ...
            'high_mean',mean(hc,'omitnan'),'high_std',std(hc,'omitnan'));
    end
    fprintf('  Budget %d%% done.\n\n',pct);
end

fid=fopen('E2_iso_summary.csv','w');
fprintf(fid,'BudgetPct,Strategy,AvgBits,HighBandAvgBits,MidDist_Mean,MidDist_Std,HighDist_Mean,HighDist_Std\n');
for b=1:N_BUD; for s=1:N_STRAT
    d=summary_data{b,s};
    fprintf(fid,'%d,%s,%.4f,%.4f,%.8f,%.8f,%.8f,%.8f\n', ...
            d.budget_pct,d.strategy,d.avg_bits,d.high_bits, ...
            d.mid_mean,d.mid_std,d.high_mean,d.high_std);
end; end
fclose(fid);
fprintf('Saved: E2_iso_summary.csv\n');

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
function d=env_dist(a,b); d=norm(a-b)/(norm(a)+eps); end

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
    fprintf('  %-12s  %6s  %6s  %8s  %10s  %12s\n', ...
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