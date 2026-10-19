function [pyc]=spinsPycnocline(p,rho_target)

params = spins_params();
dt = params.plot_interval;
x = xgrid_reader(); %x=x-max(max(x/4))-2.85;
z = zgrid_reader();
den = spins_reader('rho',p);
[Nx,Nz]=size(den);


%%% Recover pycnocline
pyc = NaN(1,Nx);

for q = 1:Nx

    rho = den(q,:);
    zz  = z(q,:);
    
    % Find where density crosses target
    ii = find((rho(1:end-1)-rho_target) .* ...
    (rho(2:end)-rho_target) <= 0, 1);
    
    if ~isempty(ii)
        
        % Linear interpolation between the two points
        pyc(q) = interp1(rho(ii:ii+1), ...
        zz(ii:ii+1), ...
        rho_target);
        
    end

end







