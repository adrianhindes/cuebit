

function [orbit] = read_CUEBIT(fname)

fid=fopen(fname,'r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;

orbit.Z       = fscanf(fid,'Z =  %g\n',1);
orbit.A       = fscanf(fid,'A =  %g\n',1);
orbit.deltat  = fscanf(fid,'delta-t = %g\n',1);
orbit.Nt      = fscanf(fid,'No. of time steps =  %d\n',1);
orbit.Ncoils  = fscanf(fid,'No. of coils = %d\n',1);
orbit.Rcoil   = fscanf(fid,'Major radius of coils =  %f\n',1);


junk          = fgets(fid);
junk          = fgets(fid);

orbit.E       = fscanf(fid,'%g %g %g',[3 inf]);
orbit.t       = orbit.deltat * (0:1:size(orbit.E,2)-1);


fclose(fid);
return;

