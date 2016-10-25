

function [orbit] = read_CUEBIT1(fname)
fname = 'full_orbit.dat'
fid=fopen(fname,'r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;



junk          = fgets(fid);

orbit.t      = fscanf(fid,'%g',1);
orbit.Zm


fclose(fid);
return;

