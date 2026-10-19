function uhat = dealias(uhat,cut)
    
    %%% Dealiase, using cutoff rules "cut"
    N = length(uhat);
    kcut = floor(cut/2*N);
    uhat(kcut+1:end-kcut) = 0;


end