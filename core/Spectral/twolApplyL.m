function Lv = twolApplyL(v,beta,Dxx,H,d,Lext,opts)
    % ---- Define productpad for odd-even
    mpoe = @(a,b) productPadWall(a,b,0,0,'odd','even',Lext);

    % ---- Construct operator L(w_t)
    if opts.Boussinesq
        Lv = v - 1/12*Dxx(mpoe(v,beta),'odd'); 
    elseif ~opts.Boussinesq
        Lv = 0;
    else
        warning('problem with opts.Boussinesq')
    end

end


