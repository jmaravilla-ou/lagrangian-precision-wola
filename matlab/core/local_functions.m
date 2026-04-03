% local_functions - shared signal processing functions for wola experiments.
%
% include these at the bottom of each experiment script, or add this
% file to your matlab path and call them directly.


function Xk = wola_analyze(x, M, win, hop)
% wola_analyze - wola analysis path, returns power spectrum per subband.
% output Xk is [M x n_frames], each entry is |fft_bin|^2.
    frame_len = 2*M;
    n_frames  = floor((length(x)-frame_len)/hop)+1;
    Xk        = zeros(M, n_frames);
    for f = 1:n_frames
        seg     = x((f-1)*hop+1:(f-1)*hop+frame_len);
        spec    = fft(seg .* win);
        Xk(:,f) = abs(spec(1:M)).^2;
    end
end


function [Xk_mag, phase] = wola_analyze_complex(x, M, win, hop)
% wola_analyze_complex - wola analysis path, returns magnitude and phase.
% phase is stored separately so synthesis can reconstruct the real signal.
    frame_len = 2*M;
    n_frames  = floor((length(x)-frame_len)/hop)+1;
    Xk_mag    = zeros(M, n_frames);
    phase     = zeros(M, n_frames);
    for f = 1:n_frames
        seg        = x((f-1)*hop+1:(f-1)*hop+frame_len);
        spec       = fft(seg .* win);
        spec_half  = spec(1:M);
        Xk_mag(:,f)= abs(spec_half);
        phase(:,f) = angle(spec_half);
    end
end


function x_out = wola_synthesize(Xk_mag, phase, M, win, hop, sig_len)
% wola_synthesize - conjugate-symmetric ifft + overlap-add reconstruction.
    frame_len = 2*M;
    n_frames  = size(Xk_mag,2);
    x_out     = zeros(sig_len + frame_len, 1);
    norm_acc  = zeros(sig_len + frame_len, 1);
    for f = 1:n_frames
        spec_half = Xk_mag(:,f) .* exp(1j * phase(:,f));
        spec_full = zeros(frame_len, 1);
        spec_full(1)       = real(spec_half(1));
        spec_full(2:M)     = spec_half(2:M);
        spec_full(M+1)     = 0;
        spec_full(M+2:end) = conj(flipud(spec_half(2:M)));
        frame     = real(ifft(spec_full));
        frame_win = frame .* win;
        idx = (f-1)*hop + (1:frame_len);
        idx = idx(idx <= length(x_out));
        x_out(idx)    = x_out(idx)    + frame_win(1:length(idx));
        norm_acc(idx) = norm_acc(idx) + win(1:length(idx)).^2;
    end
    norm_acc(norm_acc < 1e-8) = 1;
    x_out = x_out(1:sig_len) ./ norm_acc(1:sig_len);
end


function Xk_q = quantize_subbands(Xk_ref, Bv, B_max)
% quantize_subbands - midrise uniform quantizer applied per subband.
% models AND-masking the (B_max - Bk) lsbs of both mac operands.
    [M, n_frames] = size(Xk_ref);
    Xk_q = zeros(M, n_frames);
    for k = 1:M
        bk        = max(1, min(B_max, round(Bv(k))));
        scale     = 2^(bk-1);
        Xk_q(k,:) = round(Xk_ref(k,:)*scale)/scale;
    end
end


function d = envelope_distortion(env_ref, env_q)
% envelope_distortion - normalized l2 distance between two envelope vectors.
    d = norm(env_ref - env_q) / (norm(env_ref) + eps);
end


function sigma2_eff = compute_sigma2(corpus, M, WIN, HOP, FS_TARGET, S)
% compute_sigma2 - long-term subband power estimate over the full corpus.
% returns sigma2_eff = sigma2_raw .* S(k), normalized to unit peak.
    N_CLIPS    = numel(corpus);
    sigma2_sum = zeros(M,1);
    n_valid    = 0;
    fprintf('computing ltss across %d clips...\n', N_CLIPS);
    for c = 1:N_CLIPS
        try [x,fs] = audioread(corpus(c).path); catch; continue; end
        x = mean(x,2);
        if fs ~= FS_TARGET, x = resample(x,FS_TARGET,fs); end
        x = x / (max(abs(x))+eps);
        Xk = wola_analyze(x, M, WIN, HOP);
        sigma2_sum = sigma2_sum + mean(Xk,2);
        n_valid = n_valid + 1;
        if mod(c,500)==0, fprintf('  ltss: %d/%d\n',c,N_CLIPS); end
    end
    raw         = sigma2_sum/n_valid + eps;
    raw         = raw / max(raw);
    sigma2_eff  = raw .* S;
    sigma2_eff  = sigma2_eff / max(sigma2_eff);
    fprintf('ltss done.\n\n');
end
