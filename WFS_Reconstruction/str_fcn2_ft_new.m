function [D, idx] = str_fcn2_ft_new(ph, mask, delta)
% function D = str_fcn2_ft_new(ph, mask, delta)

    N = size(ph, 1);
    ph = ph .* mask;

    P = ft2(ph, delta);
    S = ft2(ph.^2, delta);
    W = ft2(mask, delta);
    delta_f = 1/(N*delta);
    w2 = ift2(W.*conj(W), delta_f);
    
    D = 2 * ift2(real(S.*conj(W)) - abs(P).^2, ...
        delta_f) ./ w2;
    
    areamask = sum(sum(mask)) * delta^2;
    maskcorr = ift2(abs(ft2(mask, delta)).^2, delta_f) / areamask;
    idx = maskcorr >= delta^2 / areamask;
    D(~idx) =  0;