import numpy as np

teff = 4948.43   
teff_err = 33

feh = -0.352
feh_err =  0.031

logg =  4.365
logg_err = 0.0651

n_samples = 1000000


feh_sys = 0.04
teff_sys = 60 
logg_sys = 0.1


teff_err_total = np.sqrt(teff_err**2 + teff_sys**2 )
feh_err_total = np.sqrt(feh_err**2 + feh_sys**2 )
logg_err_total = np.sqrt(logg_err**2 + logg_sys**2 )


print("Values after adding the systematic errors from Sousa et al. 2011")
print("Stellar effective Temperature: {0:5.0f}  +- {1:3.0f}".format(teff, teff_err_total))
print("Metallicity [Fe/H]: {0:4.2f} +- {1:4.2f}".format(feh, feh_err_total))
print("Gravity (without correction): {0:4.2f} +- {1:4.2f}".format(logg, logg_err_total))

teff_dist = np.random.normal(teff, teff_err_total, n_samples)
#logg_dist = np.random.normal(logg, logg_err_total, n_samples)


print()
print("Correction obtained from the comparison with asteroseismic gravities (Mortier et al. 2014, eq 4)")
print("This formula is valid for stars with an effective temperature between 5200 K and 7000 K")

delta_f_dist = np.random.normal(-3.89, 0.23, n_samples) * (teff_dist/10000.0) + np.random.normal(2.10, 0.14, n_samples)
med_dist = np.median(delta_f_dist)
std_dist = np.std(delta_f_dist)

logg_out = logg + med_dist
logg_out_error = np.sqrt(logg_err_total**2 + std_dist**2)

print("Correction: {0:4.2f} +/- {1:4.2f}".format(med_dist, std_dist))
print("Corrected gravity (with sys error): {0:4.2f} +/- {1:4.2f} ".format(logg_out, logg_out_error))
print()


print("Correction obtained from the comparison with photometric gravities (from transits, Mortier et al. 2014, eq 3) ")
print("This formula is valid for stars with an effective temperature between 4500 K and 7050 K")

delta_f_dist = np.random.normal(-4.57, 0.25, n_samples) * (teff/10000.0) + np.random.normal(2.59, 0.15, n_samples)
med_dist = np.median(delta_f_dist)
std_dist = np.std(delta_f_dist)

logg_out = logg + med_dist
logg_out_error = np.sqrt(logg_err_total**2 + std_dist**2)

print("Correction: {0:4.2f} +- {1:4.2f}".format(med_dist, std_dist))
print("Corrected gravity (with sys error): {0:4.2f} +/- {1:4.2f} ".format(logg_out, logg_out_error))


### add systemtics

