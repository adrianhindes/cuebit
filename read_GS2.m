

function [gs2] = read_GS2(fname)

fid=fopen(fname,'r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;

junk = fgets(fid);

gs2.header1  = fgets(fid);
gs2.header2  = fgets(fid);
gs2.header3  = fgets(fid);

junk        = fgets(fid);
junk        = fgets(fid);

gs2.r0      = fscanf(fid,'%f',1);
gs2.a       = fscanf(fid,'%f',1);
gs2.rmag    = fscanf(fid,'%f',1);
gs2.zmag    = fscanf(fid,'%f\n',1);

junk        = fgets(fid);

gs2.psimin  = fscanf(fid,'%f',1);
gs2.psedge  = fscanf(fid,'%f',1);
gs2.b0      = fscanf(fid,'%f',1);
gs2.ip      = fscanf(fid,'%f\n',1);

junk        = fgets(fid);
gs2.nfs     = fscanf(fid,'%d\n',1);

junk        = fgets(fid);
gs2.psi1D   = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.amin    = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.q       = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.f       = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.p       = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.pp      = fscanf(fid,'%f',gs2.nfs);

junk        = fgets(fid);
junk        = fgets(fid,22);
gs2.nlcfs   = fscanf(fid,'%d\n');

junk        = fgets(fid);
gs2.rlcf    = fscanf(fid,'%f',gs2.nlcfs);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.zlcf    = fscanf(fid,'%f',gs2.nlcfs);

junk        = fgets(fid);
junk        = fgets(fid);

gs2.NR      = fscanf(fid,'%d',1);
gs2.NZ      = fscanf(fid,'%d\n',1);

junk        = fgets(fid);
gs2.rgrid  = fscanf(fid,'%f',gs2.NR);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.zgrid  = fscanf(fid,'%f',gs2.NZ);

junk        = fgets(fid);
junk        = fgets(fid);
gs2.psi     = fscanf(fid,'%f',gs2.NR * gs2.NZ);
gs2.psi     = reshape(gs2.psi, gs2.NR, gs2.NZ);

[gs2.rmesh, gs2.zmesh] = meshgrid(gs2.rgrid, gs2.zgrid);
gs2.rmesh  = gs2.rmesh';
gs2.zmesh  = gs2.zmesh';

fclose(fid);
return;

