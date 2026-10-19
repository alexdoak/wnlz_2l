function fg = dealiasProduct(f,g,f_order,g_order,L)
% Function to computer delaliased product of f_{f_order} .* g_{g_order}.
% f_order is the order of the derivative for f we want in the product.
% f and g are taken from x=-L/2 to x=0, and reflected to [-L/2,L/2-dx]    
% The product is taken by double padding.
    
    % Extend data
    N0 = numel(f);
    fext = [f(1:N0);f(end-1:-1:2)];
    gext = [g(1:N0);g(end-1:-1:2)];

    N = numel(fext);
    M = 2*N;
    
    % Fourier mesh
    kx =  (2*pi/L)*[0:N/2-1 0 -N/2+1:-1]'; 

    % Fourier coefficients
    fh = (1j*kx).^f_order.*fft(fext);
    gh = (1j*kx).^g_order.*fft(gext);

    % Zero-pad in Fourier space
    fh_pad = zeros(M,1);
    gh_pad = zeros(M,1);

    % Positive frequencies including zero
    fh_pad(1:N/2) = fh(1:N/2);
    gh_pad(1:N/2) = gh(1:N/2);

    % Negative frequencies
    fh_pad(M-N/2+2:M) = fh(N/2+2:N);
    gh_pad(M-N/2+2:M) = gh(N/2+2:N);

    % Nyquist mode is omitted

    % Transform to 2x physical grid
    f_pad = (M/N)*ifft(fh_pad);
    g_pad = (M/N)*ifft(gh_pad);

    % Multiply on padded grid
    fg_pad = f_pad .* g_pad;

    % Back to Fourier space
    fgh_pad = (N/M)*fft(fg_pad);

    % Truncate back to N modes
    fgh = zeros(N,1);

    fgh(1:N/2) = fgh_pad(1:N/2);
    fgh(N/2+2:N) = fgh_pad(M-N/2+2:M);

    % Return physical-space product on original grid
    fg = real(ifft(fgh));

    fg = fg(1:N0);
end