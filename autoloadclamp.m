function [fe kv gv mv fen kvn gvn mvn] = autoloadclamp(fmin, fmax, kmin, kmax, gmin, gmax, mmin, mmax, nf, nk, ng, nm, gentypef, gentypek, gentypeg, gentypem)
% [fe kv gv mv m] = autoloadclamp(fmin, fmax, kmin, kmax, gmin, gmax, mmin, mmax, nf, nk, ng, nm, gentypef, gentypek, gentypeg, gentypem)
% gentypeX refers to linear spacing (1) or logarithhmic spacing (2) for each
% of the virtual parameters, including force


% Set up arrays 

switch gentypef
    case 1                       % linear in force
        fe = linspace(fmin,fmax,nf);
        switch gentypek
            case 1                % linear in stiffness
                kv = linspace(kmin,kmax,nk);
                switch gentypeg
                    case 1          % linear in drag
                        gv = linspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                    case 2      % logarithmic in drag
                        gv = logspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                end
            case 2              % logarithmic in stiffness
                kv = logspace(kmin,kmax,nk);
                switch gentypeg
                    case 1          % linear in drag
                        gv = linspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                    case 2      % logarithmic in drag
                        gv = logspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                end
        end
    case 2              %logarithmic in force
        fe = logspace(fmin,fmax,nf);
        switch gentypek
            case 1                % linear in stiffness
                kv = linspace(kmin,kmax,nk);
                switch gentypeg
                    case 1          % linear in drag
                        gv = linspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                    case 2      % logarithmic in drag
                        gv = logspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                end
            case 2              % logarithmic in stiffness
                kv = logspace(kmin,kmax,nk);
                switch gentypeg
                    case 1          % linear in drag
                        gv = linspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                    case 2      % logarithmic in drag
                        gv = logspace(gmin,gmax,ng);
                        switch gentypem
                            case 1      %linear in mass
                                mv = linspace(mmin,mmax,nm);
                            case 2      % logarithmic in mass
                                mv = logspace(mmin,mmax,nm);
                        end
                end
        end
end

% Concatenate arrays in order
% Order follows Fe, kv, gv, mv
r=1;
for i = 1:nf
    for j = 1:nk
        for k = 1:ng
            for l = 1:nm
                fen(r) = fe(i);
                kvn(r) = kv(j);
                gvn(r) = gv(k);
                mvn(r) = mv(l);
                r=r+1;
            end
        end
    end
end



end
