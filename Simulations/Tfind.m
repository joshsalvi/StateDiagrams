function Tfind
xfishtable = readtable('xfish1.0Noise.dat','Delimiter','\t');

time = xfishtable{3:end,{'time'}};

%Get the number of K and Fc points
Kno = xfishtable{1,{'time'}};
Fcno = xfishtable{2,{'time'}};

k_rand = zeros(Kno,Fcno);
F_rand = zeros(Kno,Fcno);
Xd = zeros(length(time),Kno,Fcno);
%File indices start at 0 but indices start at 1 in matlab
for Fcindex = 1:Fcno
    for Kindex = 1:Kno
k_rand(Kindex,Fcindex) = xfishtable{1,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
F_rand(Kindex,Fcindex) = xfishtable{2,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
Xd(:,Kindex,Fcindex) = xfishtable{3:end,{['OP' num2str(Kindex-1) num2str(Fcindex-1)]}};
    end
end

%Operating points in ascending order
Fsort = sort(F_rand(1,:));
Fgrid = Fsort(diff(Fsort) ~= 0);
Fgrid(end+1) = max(F_rand(1,:));
ksort = sort(k_rand(:,1));
kgrid = ksort(diff(ksort) ~= 0)';
kgrid(end+1) = max(k_rand(:,1));

Fgridrev = sort(Fgrid,'descend');

sizeXd = size(Xd);
Np = sizeXd(2);
Nt = sizeXd(3);

Tstartend(1:2,1:Np,1:Nt) = zeros(2,Np,Nt);

for Findex = 1:length(Fgrid)
for kindex = 1:length(kgrid)

Npt = find(k_rand == kgrid(kindex) & F_rand == Fgridrev(Findex));
if rem(Npt,Np) ~= 0;
    Npulse = rem(Npt,Np);
else
    Npulse = Np;
end  
Ntrial = (Npt - Npulse)/Np + 1;
Npt = (Ntrial-1)*Np+Npulse;

%Time for the steady state to be reached
Tstartend(1,Npulse,Ntrial) = time(1);
Tstartend(2,Npulse,Ntrial) = time(end);

end
end

%%%%%%%%%%%Save the time limits%%%%%%%%%%%
timefile = 'Tstartend.mat';
save(timefile, 'Tstartend');
display('saving...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end