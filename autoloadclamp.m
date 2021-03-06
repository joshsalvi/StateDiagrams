function [fe kv gv mv fen kvn gvn mvn] = autoloadclamp(fmin, fmax, kmin, kmax, gmin, gmax, mmin, mmax, nf, nk, ng, nm, gentypef, gentypek, gentypeg, gentypem, sorttype, sortorder)
% [fe kv gv mv m] = autoloadclamp(fmin, fmax, kmin, kmax, gmin, gmax, mmin, mmax, nf, nk, ng, nm, gentypef, gentypek, gentypeg, gentypem, sorttype, sortorder)
% gentypeX refers to linear spacing (1) or logarithhmic spacing (2) for each
% of the virtual parameters, including force
%
% Note that "sortorder" is a four-element array with unique integers 1-4.
% This defines the order of sorting, such that sortorder=[fe kv gv mv] and
% the numbers determine in which order the arrays will be sorted. For
% example sortorder=[2 1 4 3] will first sort kv, then fe, then mv, and
% finally gv.
%
% jsalvi@rockefeller.edu


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
                fkgmn(r,1) = fe(i);
                fkgmn(r,2) = kv(j);
                fkgmn(r,3) = gv(k);
                fkgmn(r,4) = mv(l);
                r=r+1;
            end
        end
    end
end

% Sorttype defines the method of iteration
% 1 = sorted; 2 = random

if sorttype == 1
    % Sort array in order specified 
    % [Fe kv gv mv] using [1 2 3 4] in any order
    fkgmn2 = sortrows(sortrows(sortrows(sortrows(fkgmn,sortorder(4)),sortorder(3)),sortorder(2)),sortorder(1));
    clear fkgmn
elseif sorttype == 2
    % Generate a random index array
    sortrand = randperm(nf*nk*ng*nm);
    for i = 1:(nf*nk*ng*nm)
        for j = 1:4
            fkgmn2(i,j) = fkgmn(sortrand(i),j);
        end
    end 
    clear fkgmn
else
    warning('No Sort Type Defined. 1=sort, 2=random');
    return;
end

% Create 1D output arrays
fen = fkgmn2(:,1);
kvn = fkgmn2(:,2);
gvn = fkgmn2(:,3);
mvn = fkgmn2(:,4);

end
