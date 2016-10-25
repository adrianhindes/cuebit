
% MJH 06/07/2011
% case with Rcoil = 99 is case without ripple
% had to edit and find/replace D+ and D- in data file with E+ and E-.

orbit1 = read_CUEBIT('MAST_loss_energies_July6_2011.dat')
orbit2 = read_CUEBIT('MAST_loss_energies1_July6_2011.dat')

plot_CUEBIT( orbit1)
plot_CUEBIT( orbit2)

[orbit] = bin_CUEBIT( orbit1)
[orbit] = bin_CUEBIT( orbit2)
