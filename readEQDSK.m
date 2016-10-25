function [eqdsk]=readEQDSK(fname)

% function reads EXPEQ and returns input quantities
% MJH 30/04/02

fid=fopen(fname,'r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;

eqdsk.comment =fscanf(fid,'%48c',1);
eqdsk.i3      =fscanf(fid,'%5d',1);
eqdsk.NRBOX   =fscanf(fid,'%5d',1);
eqdsk.NZBOX   =fscanf(fid,'%5d',1);

eqdsk.RBOXLEN =fscanf(fid,'%16e9',1);
eqdsk.ZBOXLEN =fscanf(fid,'%16e9',1);
eqdsk.R0EXP   =fscanf(fid,'%16e9',1);
eqdsk.RBOXLFT =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);

eqdsk.RAXIS   =fscanf(fid,'%16e9',1);
eqdsk.ZAXIS   =fscanf(fid,'%16e9',1);
eqdsk.PSIAXIS =fscanf(fid,'%16e9',1);
eqdsk.PSILCF  =fscanf(fid,'%16e9',1);
eqdsk.B0EXP   =fscanf(fid,'%16e9',1);

eqdsk.CURRENT =fscanf(fid,'%16e9',1);
eqdsk.PSIAXIS =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);
eqdsk.RAXIS   =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);

eqdsk.ZAXIS   =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);
junk          =fscanf(fid,'%16e9',1);

for i=1:eqdsk.NRBOX eqdsk.f(i)  =fscanf(fid,'%16e9',1); end;
for i=1:eqdsk.NRBOX eqdsk.p(i)  =fscanf(fid,'%16e9',1); end;
for i=1:eqdsk.NRBOX eqdsk.ffp(i)=fscanf(fid,'%16e9',1); end;
for i=1:eqdsk.NRBOX eqdsk.pp(i) =fscanf(fid,'%16e9',1); end;

for j=1:eqdsk.NZBOX 
 for i=1:eqdsk.NRBOX
   eqdsk.PSI(i,j) =fscanf(fid,'%16e9',1);
 end;
end;

for i=1:eqdsk.NRBOX eqdsk.q(i)  =fscanf(fid,'%16e9',1); end;

eqdsk.npbound =fscanf(fid,'%5d',1);
eqdsk.nlimiter=fscanf(fid,'%5d',1);

for i=1:eqdsk.npbound
  eqdsk.bound(i,1) =fscanf(fid,'%16e9',1);
  eqdsk.bound(i,2) =fscanf(fid,'%16e9',2); 
end;

for i=1:eqdsk.nlimiter
  eqdsk.limiter(i,1) =fscanf(fid,'%16e9',1);
  eqdsk.limiter(i,2) =fscanf(fid,'%16e9',2); 
end;

junk          =fscanf(fid,'%c',inf);

eqdsk.psn = 0:1/(eqdsk.NRBOX-1):1;
 
fclose(fid);
return;

