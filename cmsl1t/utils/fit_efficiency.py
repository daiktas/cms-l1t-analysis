from collections import OrderedDict
import math
from rootpy.plotting.func import F1
from rootpy.plotting.hist import _HistBase, Efficiency
from rootpy.plotting.graph import _GraphBase
from rootpy.ROOT import Fit, Math
from rootpy import asrootpy
import array


def fit_efficiency(efficiency, in_mean, in_sigma=10,
                   asymmetric=True, name="fit_efficiency"):
    """
    Fit a efficiency curve

    params:
    - efficiency -- a ROOT-based Efficiency plot containing the efficiency curve
                      datapoints
    - asymmetric -- use the full asymmetric fit, the cumulative dist of an
                    Exponentially Modified Gaussian (EMG), as opposed to a pure
                    Gaussian
    returns: the parameters describing the fit
    """
    efficiency = asrootpy(efficiency)
    fit_functions = []

    fit_functions.append(get_symmetric_formula())

    if asymmetric:
        fit_functions.append(get_asymmetric_formula())

    # Sometimes the mean passed in is a string, eg. when the threshold bin is "overflow" or "underflow"
    # In this case, choose a reasonable value, 50 GeV, as the starting point for the fit
    if isinstance(in_mean, str):
        in_mean = 50

    fits = []
    for i, func in enumerate(fit_functions):
        this_name = "fit_{}_{}".format(name, i)
        fitFcn = F1(func, name=this_name)
        fits.append(fitFcn)

        if i == 0:
            mu = in_mean
            sigma_inv = 1. / in_sigma
            fitFcn.SetParameters(mu, sigma_inv)
            fitFcn.SetParNames("mu", "sigma_inv")
        elif i == 1:
            mu = fits[0].GetParameter(0)
            sigma_inv = fits[0].GetParameter(1)
            lamda = 0.03  # should be within 0.04 and 0.06 it seems

            p0 = mu
            p1 = sigma_inv
            p2 = lamda
            fitFcn.SetParameters(p0, p1, p2)
            fitFcn.SetParNames("mu", "sigma_inv", "lambda_sigma")

        success = do_fit(efficiency, fitFcn)

    # Create the output parameter dictionary
    return _create_output_dict(success, fits, efficiency)


def do_fit(efficiency, fitFcn):
    # Sometimes we're given a graph, others a histogram:
    if isinstance(efficiency, _HistBase):
        opt, data = prepare_data_hist(efficiency)
    elif isinstance(efficiency, _GraphBase):
        opt, data = prepare_data_graph(efficiency)
    elif isinstance(efficiency, Efficiency):
        graph = asrootpy(efficiency.CreateGraph("e0"))
        opt, data = prepare_data_graph(graph)
    else:
        print "fit_efficiency(): Unknown object to fit of type:", type(efficiency)
        return False

    # Create the model
    fitFunction = Math.WrappedMultiTF1(fitFcn, fitFcn.GetNdim())
    fitter = Fit.Fitter()
    fitter.SetFunction(fitFunction, False)

    # Run the fit
    fitter.LikelihoodFit(data)

    # get the result
    result = fitter.Result()
    success = result.IsValid()
    return success


def prepare_data_hist(efficiency):

    # What's the range of the data to fit
    x_min = efficiency.lowerbound(0)
    x_max = efficiency.upperbound(0)
    data_range = Fit.DataRange(x_min, x_max)

    # Prepare the fit options
    opt = Fit.DataOptions()

    # Do we have y error bars?
    have_errors = False
    for a_bin in efficiency:
        if a_bin.error != 0:
            have_errors = True
            break

    # If we have errors, ignore empty bins
    if have_errors:
        opt.fUseEmpty = False
    else:
        opt.fUseEmpty = True

    data = Fit.BinData(opt, data_range)
    Fit.FillData(data, efficiency)
    return opt, data


def prepare_data_graph(efficiency):
    # What's the range of the data to fit
    x_min = efficiency.x(1) - efficiency.xerrl(1)
    x_max = efficiency.x(-2) + efficiency.xerrh(-2)
    data_range = Fit.DataRange(x_min, x_max)

    # Prepare the fit options
    opt = Fit.DataOptions()

    # Do we have y error bars?
    have_errors = False
    for i_point in range(len(efficiency)):
        if efficiency.yerrh(i_point) != 0:
            have_errors = True
            break

    # If we have errors, ignore empty bins
    if have_errors:
        opt.fUseEmpty = False
    else:
        opt.fUseEmpty = True

    data = Fit.BinData(opt, data_range)
    Fit.FillData(data, efficiency)
    return opt, data


def _create_output_dict(success, fits, input_data):
    last_fit = fits[-1]
    out_params = dict(success=success, symmetric=fits[0])
    asymmetric = len(fits) > 1
    if asymmetric:
        out_params["asymmetric"] = fits[1]
    for i in range(last_fit.GetNpar()):
        par = last_fit.GetParameter(i)
        err = last_fit.GetParError(i)
        name = last_fit.GetParName(i)
        out_params[name] = (par, err)

    if isinstance(input_data, _HistBase):
        mu_rec = out_params["mu"]
        x_min = input_data.lowerbound(0)
        x_max = input_data.upperbound(0)
        bin_widths = (x_max - x_min) / len(input_data)
        mu_rec = (mu_rec[0] + bin_widths / 2, mu_rec[1])

    sigma_inv, sigma_inv_err = out_params["sigma_inv"]
    if asymmetric:
        lambda_sigma, lambda_sigma_err = out_params["lambda_sigma"]
        out_params["sigma"] = (1. / sigma_inv, sigma_inv_err / sigma_inv ** 2)
        out_params["lambda"] = convert_to_lambda(sigma_inv, sigma_inv_err, lambda_sigma, lambda_sigma_err)
    else:
        out_params["sigma"] = (1. / sigma_inv, sigma_inv_err / sigma_inv ** 2)
    return out_params


def get_symmetric_formula():
    """
    Integral of a pure Gaussian:
    """
    return "0.5*(1+TMath::Erf((x-[0])*[1]))"


def get_asymmetric_formula():
    r"""
    Fit with an exponentially modified Gaussian (EMG)
    From Wikipedia ( https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution ) the CDF of an EMG is:
    CDF = \Phi (u,0,v)-e^{-u+v^{2}/2+\log(\Phi (u,v^{2},v))}}
    \Phi (x,\mu ,\sigma ) is the CDF of a Gaussian distribution,
    u=\lambda (x-\mu )
    v=\lambda \sigma
    \lambda (>0) := the exponential decay parameter
    \mu := the mean of the Gaussian component
    \sigma^2 (>0):= the variance of the Gaussian component

    Which simplifies to:
    func = T1 - e_m * T_2
    T1 = 0.5 * [ 1 + erf (x_prime) ]
    e_m = exp( - \lamda \sigma [ x_prime - \lamda \sigma / 2 ] )
    T2 = 0.5 * [ 1 + erf (x_prime - \lamda \sigma) ]
    x_prime = (x - \mu ) / \sigma
    so setting:
    * [0] = \mu
    * [1] = 1 / \sigma
    * [2] = \lambda \sigma
    """

    eqn = OrderedDict()
    eqn["x_prime"] = "[1]*(x - [0])"
    eqn["term_1"] = "0.5 * (1 + TMath::Erf({x_prime}))"
    eqn["exp_modif"] = "exp(-[2]*({x_prime} - [2]*0.5) )"
    eqn["term_2"] = "0.5 * (1 + TMath::Erf({x_prime} - [2]))"
    eqn["func"] = "{term_1} - {exp_modif}*{term_2}"

    for name, formula in eqn.items():
        eqn[name] = formula.format(**eqn)
    return eqn["func"]


def convert_to_lambda(p_one, p_one_err, p_two, p_two_err):
    lamda = p_two * p_one
    lamda_err = (p_two * p_one_err) ** 2
    lamda_err += (p_one * p_two_err) ** 2
    lamda_err = math.sqrt(lamda_err)
    return lamda, lamda_err
