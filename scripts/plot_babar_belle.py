import numpy as np
import argparse
from scipy.stats import multivariate_normal
import uncertainties
from uncertainties import unumpy

def setup_matplotlib(save=False):
    import matplotlib

    if save:
        # default \textwidth for a fullpage article in Latex is 16.50764 cm.
        # You can figure this out by compiling the following TeX document:
        #
        #    \documentclass{article}
        #    \usepackage{fullpage}
        #    \usepackage{layouts}
        #    \begin{document}
        #    textwidth in cm: \printinunitsof{cm}\printlen{\textwidth}
        #    \end{document}

        width = 16.50764
        width /= 2.54 # cm -> inches
        # According to this page:
        # http://www-personal.umich.edu/~jpboyd/eng403_chap2_tuftegospel.pdf,
        # Tufte suggests an aspect ratio of 1.5 - 1.6.
        height = width/1.5
        matplotlib.rcParams['figure.figsize'] = (width,height)

        matplotlib.rcParams['font.family'] = 'serif'
        matplotlib.rcParams['font.serif'] = 'computer modern roman'

        matplotlib.rcParams['text.usetex'] = True
    else:
        # on retina screens, the default plots are way too small
        # by using Qt5 and setting QT_AUTO_SCREEN_SCALE_FACTOR=1
        # Qt5 will scale everything using the dpi in ~/.Xresources
        matplotlib.use("Qt5Agg")

        # Default figure size. Currently set to my monitor width and height so that
        # things are properly formatted
        matplotlib.rcParams['figure.figsize'] = (13.78,7.48)

        # Make the defalt font bigger
        matplotlib.rcParams['font.size'] = 12

def make_cov(rho2_R1,rho2_R2,R1_R2,errors):
    """
    Returns the covariance matrix from the correlations and errors.
    """
    a = np.zeros((3,3))
    a[0, 1] = rho2_R1
    a[0, 2] = rho2_R2
    a[1, 2] = R1_R2

    corrM = a + a.T + np.identity(3)
    covM = np.atleast_2d(errors).T*corrM*errors
    return covM

def log_prob(x,mean,cov):
    return -nll(x,mean,cov)

def nll(x,mean,cov):
    return -multivariate_normal.logpdf(x,mean,cov)

if __name__ == '__main__':
    import emcee

    parser = argparse.ArgumentParser("Plot BaBar and Belle results along with fit results")
    parser.add_argument("--save", default=False, action='store_true', help="setup matplotlib to save the results")
    args = parser.parse_args()

    setup_matplotlib(save=args.save)

    import matplotlib.pyplot as plt
    import matplotlib.transforms as transforms

    names = [r'$\rho^2$', 'R1', 'R2']
    comb_par = np.array([1.122, 1.270, 0.852])
    comb_err = np.array([0.024, 0.026, 0.018])
    comb_corr = np.array([0.566, -0.824, -0.715])
    comb_cov = make_cov(*comb_corr,comb_err)
    comb_pars = uncertainties.correlated_values(comb_par,comb_cov)
    babar_par = np.array([1.191, 1.429, 0.827])
    babar_err = np.array([0.056, 0.075, 0.044])
    babar_cov = make_cov(0.71, -0.83, -0.84,babar_err)
    babar_pars = uncertainties.correlated_values(babar_par,babar_cov)
    belle_par = np.array([1.106, 1.229, 0.852])
    belle_err = np.array([0.032, 0.029, 0.022])
    belle_cov = make_cov(0.593, -0.883, -0.692,belle_err)
    belle_pars = uncertainties.correlated_values(belle_par,belle_cov)

    #samples = multivariate_normal.rvs(babar_par,babar_cov,size=100000)
    #nll_samples = np.array([nll(x,babar_par,babar_cov) for x in samples])

    #plt.hist(nll_samples,bins=100,histtype='step')
    #belle_value = nll(belle_par,babar_par,babar_cov) 
    #plt.axvline(belle_value,label="Belle Parameters",color='k',alpha=0.5,ls='--')
    #print("p-value = %.4f%%" % (100*np.count_nonzero(nll_samples > belle_value)/len(nll_samples)))
    #plt.xlabel(r"$-\log(L)$")
    #plt.title(r"$-\log(L)$ Distribution for Parameters From BaBar")
    #plt.legend()
    #plt.savefig("belle_p_value.pdf")
    #plt.savefig("belle_p_value.png")

    eigVal, eigVec = np.linalg.eig(comb_cov)
    eigSig = np.sqrt(eigVal)

    # From v16_1_noSigCuts_comb
    constrained_fit_eig = np.array([-1.08,5.31,2.46])*2
    # From v16_1_comb_CLN
    constrained_fit_eig = np.array([-1.5,3.33,0.9])*2
    constrained_fit_eig_err = np.array([0.86,0.81,0.82])*2
    constrained_fit_cov = make_cov(0.12, -0.07, -0.05,constrained_fit_eig_err)
    constrained_fit_eigs = uncertainties.correlated_values(constrained_fit_eig,constrained_fit_cov)
    # From v16_1_noSigCuts_FFunc10_comb
    unconstrained_fit_eig = np.array([-0.31,0.72,0.08])*20
    # From v16_1_FFunc10_comb
    unconstrained_fit_eig = np.array([-0.34,0.81,0.09])*20
    unconstrained_fit_eig_err = np.array([0.18,0.15,0.20])*20
    unconstrained_fit_cov = make_cov(0.45, -0.07, -0.18,unconstrained_fit_eig_err)
    unconstrained_fit_eigs = uncertainties.correlated_values(unconstrained_fit_eig,unconstrained_fit_cov)

    constrained_fit_par = comb_par + (eigSig*(constrained_fit_eigs)*eigVec).sum(axis=-1)
    unconstrained_fit_par = comb_par + (eigSig*(unconstrained_fit_eigs)*eigVec).sum(axis=-1)

    plt.figure()
    plt.subplot(2,1,1)
    offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData

    plt.errorbar(names,babar_par,yerr=babar_err,marker='o',ls='None',label="BaBar",transform=trans+offset(-20))
    plt.errorbar(names,belle_par,yerr=belle_err,marker='o',ls='None',label="Belle",transform=trans+offset(-10))
    plt.errorbar(names,comb_par,yerr=comb_err,marker='o',ls='None',label="BaBar + Belle",transform=trans+offset(0))
    plt.errorbar(names,unumpy.nominal_values(unconstrained_fit_par),yerr=unumpy.std_devs(unconstrained_fit_par),marker='o',ls='None',label="Unconstrained Fit Results",transform=trans+offset(+10))
    plt.errorbar(names,unumpy.nominal_values(constrained_fit_par),yerr=unumpy.std_devs(constrained_fit_par),marker='o',ls='None',label="Constrained Fit Results",transform=trans+offset(+20))
    plt.gca().set_ylim(top=1.8)
    plt.ylabel("Value")
    plt.title(r"$B \rightarrow D^* \mu \nu$ Form Factor Parameters")
    plt.legend(loc='upper right')
    plt.margins(x=0.1)
    plt.subplot(2,1,2)
    offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
    trans = plt.gca().transData
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.errorbar(names,unumpy.nominal_values((babar_pars-comb_par)/comb_err),yerr=unumpy.std_devs((babar_pars-comb_par)/comb_err),marker='o',ls='None',label="BaBar",transform=trans+offset(-20))
    plt.errorbar(names,unumpy.nominal_values((belle_pars-comb_par)/comb_err),yerr=unumpy.std_devs((belle_pars-comb_par)/comb_err),marker='o',ls='None',label="Belle",transform=trans+offset(-10))
    plt.errorbar(names,unumpy.nominal_values((comb_pars-comb_par)/comb_err),yerr=unumpy.std_devs((comb_pars-comb_par)/comb_err),marker='o',ls='None',label="Babar + Belle",transform=trans+offset(0))
    plt.errorbar(names,unumpy.nominal_values((unconstrained_fit_par-comb_par)/comb_err),yerr=unumpy.std_devs((unconstrained_fit_par-comb_par)/comb_err),marker='o',ls='None',label="Unconstrained",transform=trans+offset(+10))
    plt.errorbar(names,unumpy.nominal_values((constrained_fit_par-comb_par)/comb_err),yerr=unumpy.std_devs((constrained_fit_par-comb_par)/comb_err),marker='o',ls='None',label="Constrained",transform=trans+offset(+20))
    plt.axhline(y=0,ls='--',color='k',alpha=0.5)
    plt.ylabel(r"(Fit - Combined)/$\sigma$(Combined)")
    plt.legend(loc='upper right')
    plt.margins(x=0.1)
    plt.tight_layout()
    plt.savefig("babar_belle.pdf")
    plt.savefig("babar_belle.png")
    plt.show()
