function u = spectralFilter(u,parity,cut)
if nargin<3
    cut = 2/3;
end

if strcmp(parity,'even')      % even
    uext = [u; u(end-1:-1:2)];
elseif strcmp(parity,'odd')           % odd
    uext = [u; -u(end-1:-1:2)];
else
    warning('parity is either odd or even')
end

uhat = fft(uext);
uhat = dealias(uhat,cut);


uext = ifft(uhat);
u = real(uext(1:numel(u)));

end