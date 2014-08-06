
%Add missing an element to array to make it rectangular
Xd_freq_ratio(numel(Xd_freq_ratio)+1) = 0;
Xd_pow_ratio(numel(Xd_freq_ratio)+1) = 0;

%Make color map matrix
k = 1;
Fnum = 5;
for i = 1:Fnum:numel(Xd_freq_ratio)
 for j = 1:Fnum
  fmat(k,j) = Xd_freq_ratio(i+j-1);
  amat(k,j) = Xd_pow_ratio(i+j-1);
 end
 kvec(k) = k_rand(i);
 k = k+1;
end
Fvec = F_rand(1:Fnum);

%Sort matrix
[Fsvec, Findex] = sort(Fvec,'descend');
[ksvec, kindex] = sort(kvec);
fsmat = fmat(kindex,Findex)';
asmat = amat(kindex,Findex)';

%Plot matrix
imagesc(ksvec,Fsvec,fsmat);
set(gca,'YDir','normal');
colormap(jet);

figure
imagesc(ksvec,Fsvec,asmat);
set(gca,'YDir','normal');
colormap(jet);