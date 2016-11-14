fid = fopen('full_orbit.dat','r');

if fid==-1
  eqdsk=0;
  crash=1;
  return;
end;


header = fgets(fid); %read first line to remove header

orbit.data = textscan(fid,'%f %f %f %f %f %f');
orbit.t = orbit.data{1,1};
orbit.x = orbit.data{1,2};
orbit.y = orbit.data{1,3};
orbit.r = orbit.data{1,4};
orbit.z = orbit.data{1,5};
orbit.phi = orbit.data{1,6};


disp('done');

time1=10000
time2=10000


figure
plot3(orbit.x(1:time1),orbit.y(1:time1),orbit.z(1:time1))
daspect([1 1 1])
grid on

title('Orbit 3D Trace')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

figure

plot(orbit.x(1:time2),orbit.y(1:time2))
grid on

title('Orbit 2D Trace')
xlabel('x (m)')
ylabel('y (m)')
fclose(fid);

