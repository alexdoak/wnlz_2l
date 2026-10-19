function fg = productPadWall(f,g,f_order,g_order,f_parity,g_parity,Lext)
% Function to computer delaliased product of f_{f_order} .* g_{g_order}.
% f_order is the order of the derivative for f we want in the product.
% f and g are taken from x=-L/2 to x=0, and reflected to [-L/2,L/2-dx]    
% The product is taken by double padding.



    % ----- Extend f to periodic domain
    if strcmp(f_parity,'even')
        fext = [f;  f(end-1:-1:2)];
    elseif strcmp(f_parity,'odd')
        fext = [f; -f(end-1:-1:2)];
    else
        error('parity must be even or odd')
    end

    % ----- Extend g to periodic domain
    if strcmp(g_parity,'even')
        gext = [g;  g(end-1:-1:2)];
    elseif strcmp(g_parity,'odd')
        gext = [g; -g(end-1:-1:2)];
    else
        error('parity must be even or odd')
    end
    
    % ---- Tracking size of original, extension, and padded
    N0 = numel(f);
    N = numel(gext);
    M = 2*N;

    % ---- Fourier mesh
    kx =  (2*pi/Lext)*[0:N/2-1 0 -N/2+1:-1]'; 

    % ---- Fourier coefficients
    fh = (1j*kx).^f_order.*fft(fext);
    gh = (1j*kx).^g_order.*fft(gext);

    % ---- Zero-pad in Fourier space
    fh_pad = zeros(M,1);
    gh_pad = zeros(M,1);

    % ---- Positive frequencies including zero
    fh_pad(1:N/2) = fh(1:N/2);
    gh_pad(1:N/2) = gh(1:N/2);

    % ---- Negative frequencies
    fh_pad(M-N/2+2:M) = fh(N/2+2:N);
    gh_pad(M-N/2+2:M) = gh(N/2+2:N);

    % ---- Nyquist mode is omitted

    % ---- Transform to 2x physical grid
    f_pad = (M/N)*ifft(fh_pad);
    g_pad = (M/N)*ifft(gh_pad);

    % ---- Multiply on padded grid
    fg_pad = f_pad .* g_pad;

    % ---- Back to Fourier space
    fgh_pad = (N/M)*fft(fg_pad);

    % ---- Truncate back to N modes
    fgh = zeros(N,1);

    fgh(1:N/2) = fgh_pad(1:N/2);
    fgh(N/2+2:N) = fgh_pad(M-N/2+2:M);

    % ---- Return physical-space product on original grid
    fg = real(ifft(fgh));
    fg = fg(1:N0);

end

%{
function c = productPad(a,b,parity1,parity2)

    % ----- Extend a to periodic domain
    if strcmp(parity1,'even')
        aext = [a;  a(end-1:-1:2)];
    elseif strcmp(parity1,'odd')
        aext = [a; -a(end-1:-1:2)];
    else
        error('parity must be even or odd')
    end

    % ----- Extend b to periodic domain
    if strcmp(parity2,'even')
        bext = [b;  b(end-1:-1:2)];
    elseif strcmp(parity2,'odd')
        bext = [b; -b(end-1:-1:2)];
    else
        error('parity must be even or odd')
    end


    N = length(aext);

    % ----- 3/2 padding
    M = 2^nextpow2(3*N/2);
    ahat = fftshift(fft(aext));
    bhat = fftshift(fft(bext));

    pad = (M-N)/2;
    
    ahat_pad = [zeros(pad,1);
                ahat;
                zeros(pad,1)];
    bhat_pad = [zeros(pad,1);
                bhat;
                zeros(pad,1)];

    ahat_pad = ifftshift(ahat_pad);
    bhat_pad = ifftshift(bhat_pad);

    
    % ----- interpolate to finer grid
    apad = ifft(ahat_pad)*(M/N);
    bpad = ifft(bhat_pad)*(M/N);

    % ----- nonlinear product on fine grid
    cpad = apad .* bpad;

    % ----- transform back
    chat_pad = fft(cpad);

    % ----- truncate back to original spectrum
    chat_pad = fftshift(chat_pad);

    chat = chat_pad(pad+1:pad+N);
    
    chat = ifftshift(chat);
    
    chat = chat*(N/M);

    % ----- return to physical space
    cext = real(ifft(chat));

    c = cext(1:(N+2)/2);

end
%}