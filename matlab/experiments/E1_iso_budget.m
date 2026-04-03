% e1_iso_budget - perceptual and unweighted distortion vs. iso-budget.
%
% computes D_w and D_mse for uniform, linear, and proposed allocations
% at four budget levels (40/60/80/100% of the 1024-bit reference).
%
% two distortion metrics:
%   D_w   = sum_k S(k) * sigma2_raw(k) * 2^(-2*Bk)  (psychoacoustic-weighted)
%   D_mse = sum_k sigma2_raw(k) * 2^(-2*Bk)          (unweighted mse)
%
% note on sigma2_raw vs sigma2_eff:
%   sigma2_raw(k) is corpus-average subband power, used for D_w and D_mse.
%   sigma2_eff(k) = sigma2_raw(k)*S(k), normalized, fed to the lagrangian solver.
%   both are saved as csv files for use by downstream experiments.
%
% outputs:
%   E1_iso_summary.csv
%   E1_distortion_ratios.csv
%   E1_sigma2_raw.csv
%   E1_sigma2_eff.csv
%   E1_iso_bit_allocation_<pct>pct.csv  (x4)
%   E1_iso_envelope_correlation_<pct>pct.csv  (x4)
%
% run this first: all other experiments depend on E1_sigma2_raw.csv
% and E1_sigma2_eff.csv.
clear; clc; close all;

LIBRISPEECH_ROOTS = {
    '/path/to/librispeech/dev-clean', ...
    '/path/to/librispeech/dev-other'
};

% ── Parameters ────────────────────────────────────────────────────────────
M           = 64;
FS_TARGET   = 16000;
B_UNIFORM   = 16;
B_MAX   = 24;
B_MIN   = 4;
               %   B=1: sign-only quantizer, physically degenerate
               %   B=2,3: insufficient for practical audio DSP
               %   B=4: minimum for meaningful fixed-point representation
ALPHA       = 1.37;
BETA        = 0.00;
FRAME       = 2*M;
HOP         = M;
WIN         = hann(FRAME,'periodic');
C_TOTAL_REF = M * B_UNIFORM;
BUDGET_FRACS = [0.40, 0.60, 0.80, 1.00];

fc = (1:M)' * (FS_TARGET/2) / M;
S  = exp(-0.5*(log(fc/2500)/0.55).^2); S = S/max(S);
mid_bands  = find(fc>=1000 & fc<=4000);
high_bands = find(fc>4000);

strategies = {'Uniform','Linear','Quadratic'};
N_STRAT = 3; N_BUD = numel(BUDGET_FRACS);

fprintf('ISO Budget: %d bits  (M=%d x B_ref=%d)\n', C_TOTAL_REF, M, B_UNIFORM);
fprintf('Cost model: J(B) = %.2f*B^2  (beta=0)\n', ALPHA);
fprintf('Strategies: Uniform, Linear, Quadratic\n\n');

% ── Corpus ────────────────────────────────────────────────────────────────
corpus  = discover_librispeech(LIBRISPEECH_ROOTS);
N_CLIPS = numel(corpus);

% ── LTSS ──────────────────────────────────────────────────────────────────
fprintf('Computing LTSS across %d clips...\n', N_CLIPS);
sigma2_sum = zeros(M,1); n_valid = 0;
for c = 1:N_CLIPS
    try [x,fs] = audioread(corpus(c).path); catch; continue; end
    x = mean(x,2);
    if fs ~= FS_TARGET, x = resample(x,FS_TARGET,fs); end
    x = x/(max(abs(x))+eps);
    Xk = wola_analyze(x,M,WIN,HOP);
    sigma2_sum = sigma2_sum + mean(Xk,2);
    n_valid = n_valid + 1;
    if mod(c,500)==0, fprintf('  LTSS: %d/%d\n',c,N_CLIPS); end
end
fprintf('LTSS done (%d valid clips).\n\n', n_valid);

% ── Variance vectors ──────────────────────────────────────────────────────
sigma2_raw = sigma2_sum/n_valid + eps;
sigma2_raw = sigma2_raw/max(sigma2_raw);       % normalize to unit peak

sigma2_eff = sigma2_raw .* S;
sigma2_eff = sigma2_eff/max(sigma2_eff);

% Save sigma2_raw (new — enables true D_mse computation)
fid = fopen('E1_sigma2_raw.csv','w');
fprintf(fid,'BandIndex,FreqHz,SensS,sigma2_raw\n');
for k = 1:M
    fprintf(fid,'%d,%.4f,%.6f,%.12e\n',k,fc(k),S(k),sigma2_raw(k));
end
fclose(fid);
fprintf('Saved: E1_sigma2_raw.csv\n');

% Save sigma2_eff (used by E10)
fid = fopen('E1_sigma2_eff.csv','w');
fprintf(fid,'BandIndex,FreqHz,sigma2_eff\n');
for k = 1:M
    fprintf(fid,'%d,%.4f,%.12e\n',k,fc(k),sigma2_eff(k));
end
fclose(fid);
fprintf('Saved: E1_sigma2_eff.csv\n\n');

% ── Accumulators ──────────────────────────────────────────────────────────
summary_bits       = zeros(N_BUD,N_STRAT);
summary_hbits      = zeros(N_BUD,N_STRAT);
summary_mid        = zeros(N_BUD,N_STRAT);
summary_high       = zeros(N_BUD,N_STRAT);
summary_mid_std    = zeros(N_BUD,N_STRAT);
summary_high_std   = zeros(N_BUD,N_STRAT);
summary_wdist      = zeros(N_BUD,N_STRAT);  % D_w   analytical
summary_mse        = zeros(N_BUD,N_STRAT);  % D_mse analytical
summary_wdist_clip = zeros(N_BUD,N_STRAT);

% ── Main loop ─────────────────────────────────────────────────────────────
for b = 1:N_BUD
    frac    = BUDGET_FRACS(b);
    pct     = round(frac*100);
    C_total = round(frac * C_TOTAL_REF);

    fprintf('============================================================\n');
    fprintf('BUDGET: %d%%  (C_total = %d bits)\n', pct, C_total);
    fprintf('============================================================\n');

    B_u = alloc_uniform(C_total, M);
    B_l = alloc_linear(sigma2_eff, C_total, M, B_MIN, B_MAX);
    B_q = alloc_quadratic(sigma2_eff, ALPHA, BETA, C_total, M, B_MIN);
    alloc_vecs = {B_u, B_l, B_q};

    % Budget verification
    fprintf('%-12s  %6s  %6s  %6s  %8s\n','Strategy','Min','Max','Avg','Sum');
    fprintf('%s\n',repmat('-',1,45));
    for s=1:N_STRAT
        Bv=alloc_vecs{s};
        fprintf('%-12s  %6.2f  %6.2f  %6.2f  %8.1f\n',...
                strategies{s},min(Bv),max(Bv),mean(Bv),sum(Bv));
        summary_bits(b,s)  = mean(Bv);
        summary_hbits(b,s) = mean(Bv(high_bands));
    end
    fprintf('\n');

    % Allocation CSV — includes sigma2_raw for external verification
    fname = sprintf('E1_iso_bit_allocation_%dpct.csv',pct);
    fid = fopen(fname,'w');
    fprintf(fid,'BandIndex,FreqHz,SensS,sigma2_raw,Uniform,Linear,Quadratic\n');
    for k=1:M
        fprintf(fid,'%d,%.1f,%.6f,%.12e,%.4f,%.4f,%.4f\n',...
                k,fc(k),S(k),sigma2_raw(k),B_u(k),B_l(k),B_q(k));
    end
    fclose(fid);

    % ── Analytical distortion — BOTH metrics ──────────────────────────────
    % D_w   uses sigma2_raw(k) weighted by S(k)
    % D_mse uses sigma2_raw(k) without any weighting
    for s=1:N_STRAT
        Bv = alloc_vecs{s};
        summary_wdist(b,s) = sum(S          .* sigma2_raw .* 2.^(-2*Bv));
        summary_mse(b,s)   = sum(sigma2_raw              .* 2.^(-2*Bv));
    end

    fprintf('  Analytical distortion (D_w and D_mse):\n');
    u_dw  = summary_wdist(b,1);
    u_mse = summary_mse(b,1);
    fprintf('  %-12s  %16s  %16s  %12s  %12s\n', ...
            'Strategy','D_w','D_mse','U/X D_w','U/X D_mse');
    fprintf('  %s\n', repmat('-',1,72));
    for s=1:N_STRAT
        dw = summary_wdist(b,s);
        dm = summary_mse(b,s);
        if s==1
            r_dw='---'; r_dm='---';
        else
            r_dw  = sprintf('%.3fx', u_dw/dw);
            r_dm  = sprintf('%.3fx', u_mse/dm);
        end
        fprintf('  %-12s  %16.6e  %16.6e  %12s  %12s\n', ...
                strategies{s}, dw, dm, r_dw, r_dm);
    end
    fprintf('\n');

    % ── Corpus distortion ─────────────────────────────────────────────────
    corr_results = zeros(N_CLIPS,N_STRAT,3);
    fprintf('Computing corpus distortion (%d clips)...\n',N_CLIPS);
    for c = 1:N_CLIPS
        try [x,fs]=audioread(corpus(c).path); catch
            corr_results(c,:,:)=NaN; continue; end
        x=mean(x,2); if fs~=FS_TARGET, x=resample(x,FS_TARGET,fs); end
        x=x/(max(abs(x))+eps);
        Xk_ref = wola_analyze(x,M,WIN,HOP);
        for s=1:N_STRAT
            Bv=alloc_vecs{s}; Xk_q=quant_sub(Xk_ref,Bv,B_MAX);
            corr_results(c,s,1)=env_dist(rms(Xk_ref(mid_bands,:),2),...
                                         rms(Xk_q(mid_bands,:),2));
            corr_results(c,s,2)=env_dist(rms(Xk_ref(high_bands,:),2),...
                                         rms(Xk_q(high_bands,:),2));
            band_dist = zeros(M,1);
            for k=1:M
                band_dist(k) = mean((Xk_ref(k,:)-Xk_q(k,:)).^2);
            end
            corr_results(c,s,3) = sum(S .* band_dist);
        end
        if mod(c,500)==0, fprintf('  %d/%d\n',c,N_CLIPS); end
    end

    fname = sprintf('E1_iso_envelope_correlation_%dpct.csv',pct);
    fid = fopen(fname,'w');
    fprintf(fid,'Clip,Speaker,Corpus,Region,Uniform,Linear,Quadratic\n');
    reg={'Mid-Band(1-4kHz)','High-Band(>4kHz)','Weighted-Perceptual'};
    for c=1:N_CLIPS; for r=1:3
        fprintf(fid,'%s,%s,%s,%s,%.8f,%.8f,%.8f\n',...
                corpus(c).label,corpus(c).speaker,corpus(c).corpus,reg{r},...
                corr_results(c,1,r),corr_results(c,2,r),corr_results(c,3,r));
    end; end
    fclose(fid);

    for s=1:N_STRAT
        summary_mid(b,s)       = mean(corr_results(:,s,1),'omitnan');
        summary_high(b,s)      = mean(corr_results(:,s,2),'omitnan');
        summary_mid_std(b,s)   = std(corr_results(:,s,1),'omitnan');
        summary_high_std(b,s)  = std(corr_results(:,s,2),'omitnan');
        summary_wdist_clip(b,s)= mean(corr_results(:,s,3),'omitnan');
    end
    fprintf('Saved %dpct files.\n\n',pct);
end

% ── Summary CSV ───────────────────────────────────────────────────────────
fid=fopen('E1_iso_summary.csv','w');
fprintf(fid,['BudgetPct,Strategy,AvgBits,HighBandAvgBits,' ...
             'WeightedDist_Analytical,UnweightedDist_Analytical,' ...
             'WeightedDist_Corpus,MidDist_Mean,MidDist_Std,' ...
             'HighDist_Mean,HighDist_Std\n']);
for b=1:N_BUD; pct=round(BUDGET_FRACS(b)*100);
    for s=1:N_STRAT
        fprintf(fid,'%d,%s,%.4f,%.4f,%.10e,%.10e,%.10f,%.8f,%.8f,%.8f,%.8f\n',...
                pct,strategies{s},...
                summary_bits(b,s),summary_hbits(b,s),...
                summary_wdist(b,s),summary_mse(b,s),...
                summary_wdist_clip(b,s),...
                summary_mid(b,s),summary_mid_std(b,s),...
                summary_high(b,s),summary_high_std(b,s));
    end; end
fclose(fid);
fprintf('Saved: E1_iso_summary.csv\n');

% ── Clean ratio table for paper ───────────────────────────────────────────
fid=fopen('E1_distortion_ratios.csv','w');
fprintf(fid,['BudgetPct,Strategy,Dw,Dmse,' ...
             'Dw_ratio_vs_Uniform,Dmse_ratio_vs_Uniform,' ...
             'Dw_dB_vs_Uniform,Dmse_dB_vs_Uniform\n']);
for b=1:N_BUD
    pct   = round(BUDGET_FRACS(b)*100);
    u_dw  = summary_wdist(b,1);
    u_mse = summary_mse(b,1);
    for s=1:N_STRAT
        dw    = summary_wdist(b,s);
        dm    = summary_mse(b,s);
        r_dw  = u_dw/dw;
        r_dm  = u_mse/dm;
        db_dw = 10*log10(r_dw);
        db_dm = 10*log10(r_dm);
        fprintf(fid,'%d,%s,%.10e,%.10e,%.6f,%.6f,%.4f,%.4f\n',...
                pct,strategies{s},dw,dm,r_dw,r_dm,db_dw,db_dm);
    end
end
fclose(fid);
fprintf('Saved: E1_distortion_ratios.csv\n\n');

% ── Final console summary ─────────────────────────────────────────────────
fprintf('=== E1 FINAL SUMMARY ===\n\n');
fprintf('D_w   = sum_k S(k)*sigma2_raw(k)*2^(-2Bk)  [psychoacoustic-weighted]\n');
fprintf('D_mse = sum_k sigma2_raw(k)*2^(-2Bk)        [true unweighted MSE]\n\n');
for b=1:N_BUD
    pct   = round(BUDGET_FRACS(b)*100);
    u_dw  = summary_wdist(b,1);
    u_mse = summary_mse(b,1);
    fprintf('Budget %d%%:\n',pct);
    for s=1:N_STRAT
        dw=summary_wdist(b,s); dm=summary_mse(b,s);
        if s==1
            fprintf('  %-12s  D_w=%.4e  D_mse=%.4e  [baseline]\n',...
                    strategies{s},dw,dm);
        else
            fprintf('  %-12s  D_w=%.4e (%.2fx, %+.1fdB)  D_mse=%.4e (%.2fx, %+.1fdB)\n',...
                    strategies{s},dw,u_dw/dw,10*log10(u_dw/dw),...
                    dm,u_mse/dm,10*log10(u_mse/dm));
        end
    end
    fprintf('\n');
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function Xk = wola_analyze(x,M,win,hop)
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

function Bv=alloc_uniform(C_total,M)
    B_base=floor(C_total/M); rem_=C_total-B_base*M;
    Bv=B_base*ones(M,1); Bv(1:rem_)=Bv(1:rem_)+1;
end

function Bv=alloc_linear(sigma2,C_total,M,Bfloor,Bceil)
    Bv=zeros(M,1); lo=1e-20; hi=max(sigma2)*1e8;
    for i=1:600
        lam=(lo+hi)/2; raw=0.5*log2(max(eps,sigma2/lam));
        Bv=min(Bceil,max(Bfloor,raw)); cost=sum(Bv);
        if cost>C_total, lo=lam; else, hi=lam; end
        if abs(cost-C_total)/C_total<1e-11, break; end
    end
    Bv_int=max(Bfloor,round(Bv)); leftover=C_total-sum(Bv_int);
    [~,ord]=sort(sigma2,'descend'); i=1;
    while leftover~=0 && i<=M
        k=ord(i);
        if leftover>0, Bv_int(k)=Bv_int(k)+1; leftover=leftover-1;
        elseif Bv_int(k)>Bfloor, Bv_int(k)=Bv_int(k)-1; leftover=leftover+1; end
        i=i+1;
    end
    Bv=double(Bv_int);
end

function Bv=alloc_quadratic(sigma2,alpha,beta,C_bits,M,Bfloor)
    if nargin < 6, Bfloor = 4; end  % default floor matches B_MIN
    B_avg=C_bits/M; C_lut=M*(alpha*B_avg^2+beta*B_avg);
    Bv=ones(M,1); lo=1e-30; hi=1e15;
    for i=1:600
        lam=(lo+hi)/2;
        for k=1:M; Bv(k)=newton_Bk(sigma2(k),alpha,beta,lam); end
        cost=sum(alpha*Bv.^2+beta*Bv);
        if cost>C_lut, lo=lam; else, hi=lam; end
        if abs(cost-C_lut)/C_lut<1e-11, break; end
    end
    Bv=max(Bfloor,Bv);
    Bv_int=max(Bfloor,round(Bv)); leftover=C_bits-sum(Bv_int);
    [~,ord]=sort(sigma2,'descend'); i=1;
    while leftover~=0 && i<=M
        k=ord(i);
        if leftover>0, Bv_int(k)=Bv_int(k)+1; leftover=leftover-1;
        elseif Bv_int(k)>Bfloor, Bv_int(k)=Bv_int(k)-1; leftover=leftover+1; end
        i=i+1;
    end
    Bv=double(Bv_int);
end

function Bk=newton_Bk(s2,alpha,beta,lam)
    if s2<=0, Bk=1; return; end
    Bk=max(1,0.5*log2(max(eps,s2/(lam+eps))));
    for n=1:100
        f =-2*log(2)*s2*2^(-2*Bk)+lam*(2*alpha*Bk+beta);
        df= 4*log(2)^2*s2*2^(-2*Bk)+2*lam*alpha;
        dB=f/(df+eps); Bk=max(1,Bk-dB);
        if abs(dB)<1e-11, break; end
    end
end