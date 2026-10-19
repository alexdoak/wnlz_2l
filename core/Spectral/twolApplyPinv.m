function u = twolApplyPinv(f,kx,beta0)
    % ---- odd parity for w_t, use odd extension
    fext = extensionOdd(f);
    
    % ---- Fourier transform f
    fhat = fft(fext);

    % ---- Fourier sign for this linear operator
    denom = 1 + beta0*kx.^2/12;

    % ---- Apply Fourier sign
    uhat = fhat ./ denom;

    % ---- Inverse fourier transform
    uext = ifft(uhat);

    % ---- Take values in physical domain
    u = real(uext(1:length(f)));

end