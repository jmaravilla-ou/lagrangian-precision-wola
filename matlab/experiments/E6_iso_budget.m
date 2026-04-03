% e6_iso_budget - stoi evaluation across all four budget levels.
%
% evaluates stoi for uniform, linear, and proposed allocations
% across all 5567 librispeech utterances.
%
% outputs: E6_iso_stoi_<pct>pct.csv (x4), E6_iso_summary.csv
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

if ~exist('stoi','file')
    error('stoi() not found. Requires MATLAB Audio Toolbox.'); end

strategies={'Uniform','Linear','Proposed'};
N_STRAT=3; N_BUD=numel(BUDGET_FRACS);

corpus=discover_librispeech(LIBRISPEECH_ROOTS); N_CLIPS=numel(corpus);

fprintf('Computing LTSS...\n');
ss=zeros(M,1); nv=0;
for c=1:N_CLIPS
    try [x,fs]=audioread(corpus(c).path); catch; continue; end
    x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
    x=x/(max(abs(x))+eps); Xk=wola_analyze_mag(x,M,WIN,HOP);
    ss=ss+mean(Xk,2); nv=nv+1;
    if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
end
raw=ss/nv+eps; raw=raw/max(raw);
sigma2_eff=raw.*S; sigma2_eff=sigma2_eff/max(sigma2_eff);
fprintf('LTSS done.\n\n');

summary_mean=zeros(N_BUD,N_STRAT); summary_std=zeros(N_BUD,N_STRAT);

for b=1:N_BUD
    frac=BUDGET_FRACS(b); pct=round(frac*100); C_total=round(frac*C_TOTAL_REF);
    fprintf('Budget %d%% (%d bits)\n',pct,C_total);

    B_u=allocate_uniform(C_total,M);
    B_l=allocate_linear(sigma2_eff,C_total,M,B_MIN,B_MAX);
    B_q=allocate_quadratic(sigma2_eff,ALPHA,BETA,C_total,M,B_MIN,B_MAX);
    alloc_vecs={B_u,B_l,B_q};
    verify_budgets(alloc_vecs,strategies,C_total,ALPHA);

    stoi_scores=NaN(N_CLIPS,N_STRAT);
    for c=1:N_CLIPS
        try [x,fs]=audioread(corpus(c).path); catch; continue; end
        x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
        x=x/(max(abs(x))+eps); len=length(x);
        [Xk_mag,phase]=wola_analyze_complex(x,M,WIN,HOP);
        for s=1:N_STRAT
            Bv=alloc_vecs{s}; Xk_q=quant_sub(Xk_mag,Bv,B_MAX);
            x_out=wola_synthesize(Xk_q,phase,M,WIN,HOP,len);
            x_out=max(-1,min(1,x_out));
            stoi_scores(c,s)=stoi(x,x_out,FS_TARGET);
        end
        if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
    end

    fname=sprintf('E6_iso_stoi_%dpct.csv',pct);
    fid=fopen(fname,'w');
    fprintf(fid,'Clip,Speaker,Corpus,Uniform,Linear,Proposed\n');
    for c=1:N_CLIPS
        fprintf(fid,'%s,%s,%s,%.6f,%.6f,%.6f\n', ...
                corpus(c).label,corpus(c).speaker,corpus(c).corpus, ...
                stoi_scores(c,1),stoi_scores(c,2),stoi_scores(c,3));
    end
    fclose(fid);

    for s=1:N_STRAT
        v=stoi_scores(:,s); v=v(~isnan(v));
        summary_mean(b,s)=mean(v); summary_std(b,s)=std(v);
    end
    fprintf('Saved: %s\n\n',fname);
end

fid=fopen('E6_iso_summary.csv','w');
fprintf(fid,'BudgetPct,Strategy,STOI_Mean,STOI_Std\n');
for b=1:N_BUD; pct=round(BUDGET_FRACS(b)*100);
    for s=1:N_STRAT
        fprintf(fid,'%d,%s,%.6f,%.6f\n', ...
                pct,strategies{s},summary_mean(b,s),summary_std(b,s));
    end
end
fclose(fid);
fprintf('Saved: E6_iso_summary.csv\n');

% ── LOCAL FUNCTIONS ───────────────────────────────────────────────────────
function Xk=wola_analyze_mag(x,M,win,hop)
    fl=2*M; nf=floor((length(x)-fl)/hop)+1; Xk=zeros(M,nf);
    for f=1:nf
        seg=x((f-1)*hop+1:(f-1)*hop+fl);
        sp=fft(seg.*win); Xk(:,f)=abs(sp(1:M)).^2;
    end
end
function [mag,ph]=wola_analyze_complex(x,M,win,hop)
    fl=2*M; nf=floor((length(x)-fl)/hop)+1;
    mag=zeros(M,nf); ph=zeros(M,nf);
    for f=1:nf
        seg=x((f-1)*hop+1:(f-1)*hop+fl);
        sp=fft(seg.*win); sh=sp(1:M);
        mag(:,f)=abs(sh); ph(:,f)=angle(sh);
    end
end
function x_out=wola_synthesize(mag,ph,M,win,hop,slen)
    fl=2*M; nf=size(mag,2);
    x_out=zeros(slen+fl,1); na=zeros(slen+fl,1);
    for f=1:nf
        sh=mag(:,f).*exp(1j*ph(:,f));
        sf=zeros(fl,1);
        sf(1)=real(sh(1)); sf(2:M)=sh(2:M); sf(M+1)=0;
        sf(M+2:end)=conj(flipud(sh(2:M)));
        fr=real(ifft(sf)).*win;
        idx=(f-1)*hop+(1:fl); idx=idx(idx<=length(x_out));
        x_out(idx)=x_out(idx)+fr(1:length(idx));
        na(idx)=na(idx)+win(1:length(idx)).^2;
    end
    na(na<1e-8)=1; x_out=x_out(1:slen)./na(1:slen);
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