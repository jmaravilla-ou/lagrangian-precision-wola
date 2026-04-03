% e10_iso_distortion - iso-distortion bit budget analysis.
%
% binary search for the minimum bit budget each strategy needs to match
% uniform's D_w at each quality target level.
%
% note on sigma2_raw vs sigma2_eff:
%   D_w evaluation uses sigma2_raw (true subband power).
%   the lagrangian solver receives sigma2_eff = sigma2_raw * S(k).
%   both are loaded from the E1 csv outputs.
%
% outputs: E10_iso_distortion_summary.csv, E10_iso_distortion_detail.csv
%
% run after E1.
clear; clc; close all;

M=64; FS_TARGET=16000; B_UNIFORM=16; B_MAX=24; B_MIN=4;
ALPHA=1.37; BETA=0.00;
FRAME=2*M; HOP=M; WIN=hann(FRAME,'periodic');
C_TOTAL_REF=M*B_UNIFORM;
fc=(1:M)'*(FS_TARGET/2)/M;
S=exp(-0.5*(log(fc/2500)/0.55).^2); S=S/max(S);

% ── Load D_w targets from E1 summary ─────────────────────────────────────
% D_w targets from E1 summary
if ~isfile('E1_iso_summary.csv')
    error('E1_iso_summary.csv not found. Run E1 first.'); end
T=readtable('E1_iso_summary.csv','TextType','string','VariableNamingRule','preserve');
BUDGET_PCTS     = [40,60,80,100];
UNIFORM_BUDGETS = round([0.40,0.60,0.80,1.00]*C_TOTAL_REF);
BUDGET_LABELS   = {'40%','60%','80%','100%'};
DW_TARGETS = zeros(4,1);
for b=1:4
    row=T(T.BudgetPct==BUDGET_PCTS(b) & strcmp(T.Strategy,'Uniform'),:);
    if isempty(row), error('Missing Uniform row for %d%%',BUDGET_PCTS(b)); end
    DW_TARGETS(b)=row.WeightedDist_Analytical;
end
fprintf('D_w targets (from E1, computed with sigma2_raw):\n');
for b=1:4, fprintf('  %s: %.6e\n',BUDGET_LABELS{b},DW_TARGETS(b)); end
fprintf('\n');

% ── Load sigma2_raw (for D_w evaluation) ─────────────────────────────────
% D_w = sum(S(k) * sigma2_raw(k) * 2^(-2Bk))
% sigma2_raw is the raw corpus-average subband power before S(k) weighting.
if ~isfile('E1_sigma2_raw.csv')
    error(['E1_sigma2_raw.csv not found. Run E1 first.\n' ...
           '(The updated E1 saves both sigma2_raw and sigma2_eff.)']);
end
Tr = readtable('E1_sigma2_raw.csv','VariableNamingRule','preserve');
sigma2_raw = Tr.sigma2_raw;
fprintf('Loaded sigma2_raw (%d bands) from E1_sigma2_raw.csv.\n', M);

% ── Load sigma2_eff (for allocator input) ────────────────────────────────
% sigma2_eff = sigma2_raw * S(k) / max(sigma2_raw * S(k))
% This is what the Lagrangian solver receives.
if isfile('E1_sigma2_eff.csv')
    Te = readtable('E1_sigma2_eff.csv','VariableNamingRule','preserve');
    sigma2_eff = Te.sigma2_eff;
    fprintf('Loaded sigma2_eff (%d bands) from E1_sigma2_eff.csv.\n\n', M);
else
    fprintf('E1_sigma2_eff.csv not found — deriving from sigma2_raw.\n');
    sigma2_eff = sigma2_raw .* S;
    sigma2_eff = sigma2_eff / max(sigma2_eff);
    fprintf('\n');
end

% ── Verification: recompute D_w for Uniform and compare to E1 targets ────
fprintf('=== VERIFICATION (should match E1 targets within rounding) ===\n');
for t=1:4
    Bv_u = allocate_uniform(UNIFORM_BUDGETS(t), M);
    dw_rc = compute_dw(sigma2_raw, S, Bv_u);   % uses sigma2_raw ✓
    rel = abs(dw_rc - DW_TARGETS(t)) / DW_TARGETS(t);
    flag = ''; if rel > 0.001, flag = '  <-- WARNING >0.1%'; end
    fprintf('  %s: target=%.6e  recomp=%.6e  err=%.4f%%%s\n', ...
            BUDGET_LABELS{t}, DW_TARGETS(t), dw_rc, rel*100, flag);
end
fprintf('\n');

% ── Binary search for minimum budget per strategy ────────────────────────
strat_names={'Uniform','Linear','Proposed'};
N_STRAT=3; N_TARGETS=4;
min_budget  = zeros(N_TARGETS, N_STRAT);
achieved_dw = zeros(N_TARGETS, N_STRAT);

fprintf('=== BINARY SEARCH ===\n\n');
for t=1:N_TARGETS
    target   = DW_TARGETS(t);
    C_uniform= UNIFORM_BUDGETS(t);
    fprintf('Target %s  (D_w <= %.6e, Uniform needs %d bits)\n', ...
            BUDGET_LABELS{t}, target, C_uniform);

    for s=1:N_STRAT
        % Check whether the strategy can reach the target at C_uniform
        Bv_check = run_alloc(s, sigma2_eff, C_uniform, M, ALPHA, BETA, B_MIN, B_MAX);
        if compute_dw(sigma2_raw, S, Bv_check) > target
            % Need more bits — expand upper bound
            C_hi = C_TOTAL_REF * 2;
            if compute_dw(sigma2_raw, S, ...
                    run_alloc(s,sigma2_eff,C_hi,M,ALPHA,BETA,B_MIN,B_MAX)) > target
                fprintf('  %-10s: cannot reach target\n', strat_names{s});
                min_budget(t,s) = NaN; achieved_dw(t,s) = NaN; continue;
            end
        else
            C_hi = C_uniform;
        end

        % Binary search: find smallest C s.t. D_w(strategy@C) <= target
        lo = M*B_MIN;  hi = C_hi;
        while hi - lo > 1
            mid = floor((lo+hi)/2);
            Bv_mid = run_alloc(s, sigma2_eff, mid, M, ALPHA, BETA, B_MIN, B_MAX);
            if compute_dw(sigma2_raw, S, Bv_mid) <= target
                hi = mid;
            else
                lo = mid;
            end
        end
        Bv_opt  = run_alloc(s, sigma2_eff, hi, M, ALPHA, BETA, B_MIN, B_MAX);
        dw_opt  = compute_dw(sigma2_raw, S, Bv_opt);
        min_budget(t,s)  = hi;
        achieved_dw(t,s) = dw_opt;
        sav  = C_uniform - hi;
        spct = sav / C_uniform * 100;
        fprintf('  %-10s: %4d bits  saving=%+d(%+.1f%%)  D_w=%.4e  OK=%d\n', ...
                strat_names{s}, hi, sav, spct, dw_opt, dw_opt<=target);
    end
    fprintf('\n');
end

% ── Output CSVs ──────────────────────────────────────────────────────────
fid = fopen('E10_iso_distortion_summary.csv','w');
fprintf(fid,'QualityTarget,UniformBudget_bits,');
for s=1:N_STRAT
    fprintf(fid,'%s_MinBudget,%s_Saving_bits,%s_Saving_pct,%s_AchievedDw,', ...
            strat_names{s},strat_names{s},strat_names{s},strat_names{s});
end
fprintf(fid,'Dw_Target\n');
for t=1:N_TARGETS
    fprintf(fid,'%s,%d,',BUDGET_LABELS{t},UNIFORM_BUDGETS(t));
    for s=1:N_STRAT
        mb=min_budget(t,s); adw=achieved_dw(t,s);
        sav=UNIFORM_BUDGETS(t)-mb; spct=sav/UNIFORM_BUDGETS(t)*100;
        fprintf(fid,'%d,%d,%.4f,%.6e,',mb,sav,spct,adw);
    end
    fprintf(fid,'%.6e\n',DW_TARGETS(t));
end
fclose(fid); fprintf('Saved: E10_iso_distortion_summary.csv\n');

fid=fopen('E10_iso_distortion_detail.csv','w');
fprintf(fid,'QualityTarget,Strategy,BandIndex,FreqHz,SensS,Bk\n');
for t=1:N_TARGETS; for s=1:N_STRAT
    if isnan(min_budget(t,s)), continue; end
    Bv=run_alloc(s,sigma2_eff,min_budget(t,s),M,ALPHA,BETA,B_MIN,B_MAX);
    for k=1:M
        fprintf(fid,'%s,%s,%d,%.1f,%.6f,%d\n', ...
                BUDGET_LABELS{t},strat_names{s},k,fc(k),S(k),Bv(k));
    end
end; end
fclose(fid); fprintf('Saved: E10_iso_distortion_detail.csv\n');

% ── LOCAL FUNCTIONS ───────────────────────────────────────────────────────
function dw = compute_dw(sigma2_raw, S, Bv)
    % CORRECT: uses sigma2_raw and S(k) separately, not sigma2_eff.
    % D_w = sum_k S(k) * sigma2_raw(k) * 2^(-2*Bk)
    dw = sum(S .* sigma2_raw .* 2.^(-2*Bv));
end

function Bv = run_alloc(s, sigma2_eff, C_total, M, alpha, beta, B_floor, B_ceil)
    switch s
        case 1, Bv = allocate_uniform(C_total, M);
        case 2, Bv = allocate_linear(sigma2_eff, C_total, M, B_floor, B_ceil);
        case 3, Bv = allocate_quadratic(sigma2_eff, alpha, beta, C_total, M, B_floor, B_ceil);
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