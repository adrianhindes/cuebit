fid = fopen('full_orbit.dat','r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;


header = fgets(fid); %read first line to remove header

orbit.data = textscan(fid,'%f %f %f %f');
orbit.t = orbit.data{1,1};
orbit.r = orbit.data{1,2};
orbit.z = orbit.data{1,3};
orbit.phi = orbit.data{1,4};


disp('done');

range=20000;

figure
plot(orbit.r(1:range),orbit.z(1:range))

title('Orbit Trace')
xlabel('Radius (m)')
ylabel('Z (m)')


fclose(fid);

