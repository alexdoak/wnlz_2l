function [hx,hxx,hxxx] = SpatialDifferentiation(h,N,dx,Dx,Dxx,Dxxx,kx,differentiation)

switch differentiation
    case 'spectral'
        % Reflect
        h = [h(1:N/2+1);h(end-1:-1:2)];

        % Differentiate in Fourier Space
        h_hat=fft(h);
        hx = real(ifft(1j*kx.*h_hat));
        hxx = real(ifft(-kx.^2.*h_hat));
        hxxx = real(ifft(-1j*kx.^3.*h_hat));

        % Take only half a wavelength (opposite of reflect)
        h = h(1:N/2+1);
        hx = hx(1:N/2+1);
        hxx = hxx(1:N/2+1);
        hxxx = hxxx(1:N/2+1);

    case 'findif'
        % Finite difference formula for h
        hx = Dx*h;
        hxx = Dxx*h;
        hxxx = Dxxx*h;
end


