import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc # latex in matplotlib
rc('text',usetex=True)

def log_Cubic_interpolate(xx,yy):
    from scipy.interpolate import CubicSpline
    logx=np.log10(xx)
    logy=np.log10(yy)
    Cubic_interp=CubicSpline(logx,logy)
    log_Cubic_interp=lambda zz: np.power(10, Cubic_interp(np.log10(zz)))
    return log_Cubic_interp

mass_len=50
mass_start=1e-5
mass_end=5e3
dm_mass_set=np.geomspace(mass_start,mass_end,mass_len)

sig_pann_cvl=np.load('../data/sigma_data/CVL_fsky1_lmin1_lmaxells.npy')
sig_pann_ali_nls_nofgres1000=np.load('../data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax1000.npy')
sig_pann_ali_nls_nofgres3000=np.load('../data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax3000.npy')
sig_pann_ali_nls_nofgres4000=np.load('../data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax4000.npy')

sig_pann_pico_nls_nofgres1000=np.load('../data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax1000.npy')
sig_pann_pico_nls_nofgres3000=np.load('../data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax3000.npy')
sig_pann_pico_nls_nofgres4000=np.load('../data/sigma_data/pico_noise_nofgres_fsky37_lmin10_lmax4000.npy')

sig_pann_ali_nls_nofgres620=np.load('../data/sigma_data/ali_noise_nofgres_fsky40_lmin30_lmax620.npy')
sig_pann_ali_nls_fgres620=np.load('../data/sigma_data/ali_noise_fgres_fsky37_lmin30_lmax620.npy')

xs=np.geomspace(mass_start,mass_end,1000)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_cvl[:,3])
plt.loglog(xs,cs(xs),color='dimgrey',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_ali_nls_nofgres1000[:,3])
plt.loglog(xs,cs(xs),color='orange',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_ali_nls_nofgres3000[:,3])
plt.loglog(xs,cs(xs),color='orange',linestyle=':',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_ali_nls_nofgres4000[:,3])
plt.loglog(xs,cs(xs),color='orange',linestyle='--',linewidth=1)

cs=log_Cubic_interpolate(dm_mass_set, sig_pann_pico_nls_nofgres1000[:,3])
plt.loglog(xs,cs(xs),color='blue',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_pico_nls_nofgres3000[:,3])
plt.loglog(xs,cs(xs),color='blue',linestyle=':',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_pico_nls_nofgres4000[:,3])
plt.loglog(xs,cs(xs),color='blue',linestyle='--',linewidth=1)

cs=log_Cubic_interpolate(dm_mass_set, sig_pann_ali_nls_nofgres620[:,3])
plt.loglog(xs,cs(xs),color='pink',linestyle='-',linewidth=1)
cs=log_Cubic_interpolate(dm_mass_set, sig_pann_ali_nls_fgres620[:,3])
plt.loglog(xs,cs(xs),color='purple',linestyle='-',linewidth=1)

plt.legend(['CVL fsky=1 l:1-4000','ali noise fsky=0.4 l:30-1000','ali noise fsky=0.4 l:10-3000','ali noise fsky=0.4 l:10-4000','pico noise fsky=0.7 l:10-1000','pico noise fsky=0.7 l:10-3000','pico noise fsky=0.7 l:10-4000','ali noise fsky=0.4 l:30-620','ali noise +fgres fsky=0.4 l:30-620'],loc='upper right')
# plt.legend(['Ali with fg residual fsky=0.5','ali without fg residual fsky=0.5','pico with fg residual fsky=1','pico without fg residual fsky=1','CVL fsky=1'],loc='upper right')

# plt.title('ali l from 30 to 620, pico l from 10 to 1000')
plt.xlabel(r'$m_\chi\left[GeV\right]$')
plt.ylabel(r'$\left\langle \sigma v\right\rangle/ m_\chi \left[cm^3s^{-1}GeV^{-1}\right] $')
plt.ylim((1e-28,1e-26))
plt.xlim((1e-5,5e3))
plt.title(r'constraints derived by fisher matrix')
plt.savefig('./figure/constraints_diff_mass_without_planck/1.png',dpi=300)
# plt.show()


