
% MJH 06/07/2011
% case with Rcoil = 99 is case without ripple
% had to edit and find/replace D+ and D- in data file with E+ and E-.

orbit1 = read_CUEBIT('orbit.dat')
orbit2 = read_CUEBIT('full_orbit.dat')

plot_CUEBIT( orbit1)
plot_CUEBIT( orbit2)

[orbit] = bin_CUEBIT( orbit1)
[orbit] = bin_CUEBIT( orbit2)
