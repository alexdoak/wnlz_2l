function ux = DxSpectral(u,kx,parity,order,cut)

if strcmp(parity,'even')      % even
    uext = [u; u(end-1:-1:2)];
elseif strcmp(parity,'odd')           % odd
    uext = [u; -u(end-1:-1:2)];
else
    warning('parity is either odd or even')
end

uhat = fft(uext);

if nargin>4
    uhat = dealias(uhat,cut);
end

uxext = ifft((1i*kx).^order .* uhat);

Next = numel(uext);
ux = real(uxext(1:numel(u)));

end