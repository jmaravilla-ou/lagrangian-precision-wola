% iso_budget_config - shared configuration for all iso-budget experiments.
%
% cost model: J(B) = alpha * B^2, alpha = 1.37, R^2 = 0.9998
% measured on Artix-7 FPGAs with DSP48 inference suppressed.
%
% returns a config struct used by all experiment scripts.

function cfg = iso_budget_config()

    cfg.M         = 64;
    cfg.FS_TARGET = 16000;
    cfg.B_REF     = 16;
    cfg.B_MAX     = 24;
    cfg.B_MIN     = 4;

    cfg.ALPHA = 1.37;
    cfg.BETA  = 0.00;

    cfg.C_TOTAL_REF = cfg.M * cfg.B_REF;

    cfg.BUDGET_FRACTIONS = [0.40, 0.60, 0.80, 1.00];

    fc      = (1:cfg.M)' * (cfg.FS_TARGET / 2) / cfg.M;
    S       = exp(-0.5 * (log(fc / 2500) / 0.55).^2);
    cfg.S   = S / max(S);
    cfg.fc  = fc;

    cfg.mid_bands  = find(fc >= 1000 & fc <= 4000);
    cfg.high_bands = find(fc >  4000);
    cfg.low_bands  = find(fc <  1000);

    fprintf('config: J(B) = %.2f*B^2, C_ref = %d bits, B_min = %d, B_max = %d\n\n', ...
            cfg.ALPHA, cfg.C_TOTAL_REF, cfg.B_MIN, cfg.B_MAX);
end
