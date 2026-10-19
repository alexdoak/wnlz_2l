function c = productPad(a,b,order1,order2,parity1,parity2,L)
% Function to computer delaliased product of f_{f_order} .* g_{g_order}.
% f_order is the order of the derivative for f we want in the product.
% f and g are taken from x=-L/2 to x=0, and reflected to [-L/2,L/2-dx]    
% The product is taken by double padding.


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