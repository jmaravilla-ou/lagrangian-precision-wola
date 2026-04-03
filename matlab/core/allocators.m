% allocators - bit allocation strategies for the wola filterbank.
%
% cost model: J(B) = alpha * B^2, alpha = 1.37
% measured on artix-7 fpgas with lut-only partial-product logic.
%
% strategies:
%   allocate_uniform    - floor(C_total/M) bits per band
%   allocate_linear     - water-filling, sum(Bk) = C_total
%   allocate_quadratic  - psychoacoustic lagrangian with hardware cost model
%
% all allocators enforce B_floor <= Bk <= B_ceil.
% use B_floor = 4, B_ceil = 24 to match the paper.


function Bv = allocate_uniform(C_total, M)
% allocate_uniform - equal bits per band, remainder to lowest-indexed bands.
    B_base   = floor(C_total / M);
    rem_bits = C_total - B_base * M;
    Bv       = B_base * ones(M, 1);
    Bv(1:rem_bits) = Bv(1:rem_bits) + 1;
end


function Bv = allocate_linear(sigma2, C_total, M, B_floor, B_ceil)
% allocate_linear - water-filling over the psychoacoustic objective.
% sigma2 is sigma2_eff (S(k)-weighted subband power) from the ltss step.
    Bv = zeros(M, 1);
    lo = 1e-20;  hi = max(sigma2) * 1e8;
    for iter = 1:600
        lam  = (lo + hi) / 2;
        raw  = 0.5 * log2(max(eps, sigma2 / lam));
        Bv   = min(B_ceil, max(B_floor, raw));
        cost = sum(Bv);
        if cost > C_total,  lo = lam;  else,  hi = lam;  end
        if abs(cost - C_total) / C_total < 1e-11,  break;  end
    end
    Bv_int   = min(B_ceil, max(B_floor, round(Bv)));
    leftover = C_total - sum(Bv_int);
    [~, ord] = sort(sigma2, 'descend');
    i = 1;
    while leftover ~= 0 && i <= M
        k = ord(i);
        if leftover > 0 && Bv_int(k) < B_ceil
            Bv_int(k) = Bv_int(k) + 1;
            leftover  = leftover - 1;
        elseif leftover < 0 && Bv_int(k) > B_floor
            Bv_int(k) = Bv_int(k) - 1;
            leftover  = leftover + 1;
        end
        i = i + 1;
    end
    Bv = double(Bv_int);
end


function Bv = allocate_quadratic(sigma2, alpha, beta, C_bits, M, B_floor, B_ceil)
% allocate_quadratic - lagrangian allocation with hardware cost constraint.
%
% minimizes D_w = sum(S(k) * sigma2(k) * 2^(-2Bk)) subject to
% sum(alpha * Bk^2) <= C_lut, where C_lut is the lut cost of the uniform
% allocation at the same bit count.
%
% continuous solution via newton/bisection on the kkt condition,
% then integer rounding with exact bit-count enforcement.
    if nargin < 6,  B_floor = 4;   end
    if nargin < 7,  B_ceil  = 24;  end

    B_avg = C_bits / M;
    C_lut = M * (alpha * B_avg^2 + beta * B_avg);

    Bv = B_floor * ones(M, 1);
    lo = 1e-30;  hi = 1e15;
    for iter = 1:600
        lam = (lo + hi) / 2;
        for k = 1:M
            Bv(k) = solve_Bk_newton(sigma2(k), alpha, beta, lam, B_floor);
        end
        cost = sum(alpha * Bv.^2 + beta * Bv);
        if cost > C_lut,  lo = lam;  else,  hi = lam;  end
        if abs(cost - C_lut) / C_lut < 1e-11,  break;  end
    end

    Bv = min(B_ceil, max(B_floor, Bv));

    Bv_int   = min(B_ceil, max(B_floor, round(Bv)));
    leftover = C_bits - sum(Bv_int);
    [~, ord] = sort(sigma2, 'descend');
    i = 1;
    while leftover ~= 0 && i <= M
        k = ord(i);
        if leftover > 0 && Bv_int(k) < B_ceil
            Bv_int(k) = Bv_int(k) + 1;
            leftover  = leftover - 1;
        elseif leftover < 0 && Bv_int(k) > B_floor
            Bv_int(k) = Bv_int(k) - 1;
            leftover  = leftover + 1;
        end
        i = i + 1;
    end
    Bv = double(Bv_int);
end


function Bk = solve_Bk_newton(s2, alpha, beta, lam, B_floor)
% solve_Bk_newton - newton solver for the per-band kkt condition.
%
% solves: -2*ln2*s2*2^(-2Bk) + lam*(2*alpha*Bk + beta) = 0
% initial guess is the water-filling solution.
    if nargin < 5,  B_floor = 4;  end
    if s2 <= 0,  Bk = B_floor;  return;  end
    Bk = max(B_floor, 0.5 * log2(max(eps, s2 / (lam + eps))));
    for n = 1:100
        f  = -2*log(2)*s2*2^(-2*Bk) + lam*(2*alpha*Bk + beta);
        df =  4*log(2)^2*s2*2^(-2*Bk) + 2*lam*alpha;
        dB = f / (df + eps);
        Bk = max(B_floor, Bk - dB);
        if abs(dB) < 1e-11,  break;  end
    end
end


function verify_budgets(alloc_vecs, strategies, C_total, alpha)
% verify_budgets - print allocation summary and flag any constraint violations.
    if nargin < 4,  alpha = 1.37;  end
    C_lut_uniform = alpha * sum(alloc_vecs{1}.^2);
    fprintf('%-20s  %8s  %8s  %8s  %10s  %12s  %12s  %s\n', ...
            'strategy','min(Bk)','max(Bk)','avg(Bk)', ...
            'sum(Bk)','sum(J(Bk))','C_lut_unif','ok?');
    fprintf('%s\n', repmat('-', 1, 90));
    for s = 1:numel(strategies)
        Bv     = alloc_vecs{s};
        sum_Bk = sum(Bv);
        sum_J  = alpha * sum(Bv.^2);
        ok_str = '';
        if sum_Bk ~= C_total,  ok_str = [ok_str ' sum_mismatch'];  end
        if any(Bv < 4),        ok_str = [ok_str ' floor_violated'];  end
        if any(Bv > 24),       ok_str = [ok_str ' ceil_violated'];   end
        if isempty(ok_str),    ok_str = 'ok';  end
        fprintf('%-20s  %8.0f  %8.0f  %8.2f  %10.0f  %12.2f  %12.2f  %s\n', ...
                strategies{s}, min(Bv), max(Bv), mean(Bv), ...
                sum_Bk, sum_J, C_lut_uniform, ok_str);
    end
    fprintf('\n');
end
