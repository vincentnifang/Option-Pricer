__author__ = 'vincent'

from scipy import stats
import numpy, math

STANDARD = 'Standard'
GEO_MEAN = 'Geometric mean Asian'
GEO_MEAN_STRIKE = 'Geometric mean Asian with adjusted strike'


STANDARD = 'Standard'
GEO_MEAN = 'Geometric mean Asian'
GEO_MEAN_STRIKE = 'Geometric mean Asian with adjusted strike'


def bs(s, k, t, v, r, option_type):
    """ Black-Scholes model.
    s: Stock price
    k: strike price
    t: maturity day
    v: volatility
    r: risk-free rate
    option_type: call/put
    """

    s = float(s)
    k = float(k)
    t = float(t)
    v = float(v)
    r = float(r)

    d1 = (math.log(s / k) + r * t) / (v * math.sqrt(t)) + 0.5 * v * math.sqrt(t)
    d2 = d1 - v * math.sqrt(t)

    if option_type == 'call':
        return s * stats.norm.cdf(d1) - k * math.exp(-r * t) * stats.norm.cdf(d2)
    elif option_type == 'put':
        return k * math.exp(-r * t) * stats.norm.cdf(-d2) - s * stats.norm.cdf(-d1)
    return None


def _get_geometric_sigma(V, N):
    return V * math.sqrt((N + 1) * (2 * N + 1) / (6 * N * N))


def _get_geometric_miu(R, V, N, sigma):
    return (R - 0.5 * V * V) * (N + 1) / (2 * N) + 0.5 * sigma * sigma


def geometric_asian_option(K, T, R, V, S0, N, option_type):
    sigma = _get_geometric_sigma(V, N)
    miu = _get_geometric_miu(R, V, N, sigma)
    d1 = (math.log(S0 / K) + (miu + 0.5 * sigma * sigma) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if option_type == 'call':
        return math.exp(-R * T) * (S0 * math.exp(miu * T) * stats.norm.cdf(d1) - K * stats.norm.cdf(d2))
    elif option_type == 'put':
        return math.exp(-R * T) * (K * stats.norm.cdf(-d2) - S0 * math.exp(miu * T) * stats.norm.cdf(-d1))
    return None


def geometric_basket_option(S1, S2, V1, V2, R, T, K, rou, option_type):
    bg_0 = math.sqrt(S1 * S2)
    sigma_bg = math.sqrt(V1 * V1 + V1 * V2 * rou + V2 * V1 * rou + V2 * V2) / 2
    miu_bg = R - 0.5 * (V1 * V1 + V2 * V2) / 2 + 0.5 * sigma_bg * sigma_bg

    d1 = (math.log(bg_0 / K) + (miu_bg + 0.5 * sigma_bg * sigma_bg) * T) / (sigma_bg * math.sqrt(T))
    d2 = d1 - sigma_bg * math.sqrt(T)

    if option_type == 'call':
        return math.exp(-R * T) * (bg_0 * math.exp(miu_bg * T) * stats.norm.cdf(d1) - K * stats.norm.cdf(d2))
    elif option_type == 'put':
        return math.exp(-R * T) * (K * stats.norm.cdf(-d2) - bg_0 * math.exp(miu_bg * T) * stats.norm.cdf(-d1))
    return None


def arithmetic_asian_option(K, geo_K, T, R, V, S0, N, option_type, path_num=10000, control_variate='Standard'):
    if control_variate == STANDARD:
        dt = T / N
        sigma = V
        drift = math.exp((R - 0.5 * sigma * sigma) * dt)
        # arith_payoff = []
        arith_payoff = numpy.empty(path_num)
        sigma_sqrt = sigma * math.sqrt(dt)
        exp_RT = math.exp(-R * T)
        for i in xrange(path_num):
            # spath = []
            spath = numpy.empty(N)
            former = S0
            for j in xrange(int(N)):
                # growth_factor = drift * math.exp(sigma * math.sqrt(dt) * numpy.random.normal(0, 1))
                growth_factor = drift * math.exp(sigma_sqrt * numpy.random.normal(0, 1))
                former = former * growth_factor
                # spath.append(former)
                spath[j] = former

            arith_mean = numpy.mean(spath)

            if option_type == 'call':
                # arith_payoff_call = math.exp(-R * T) * max(arith_mean - K, 0)
                arith_payoff_call = exp_RT * max(arith_mean - K, 0)
                # arith_payoff.append(arith_payoff_call)
                arith_payoff[i] = arith_payoff_call
            elif option_type == 'put':
                # arith_payoff_put = math.exp(-R * T) * max(K - arith_mean, 0)
                arith_payoff_put = exp_RT * max(K - arith_mean, 0)
                # arith_payoff.append(arith_payoff_put)
                arith_payoff[i] = arith_payoff_put


        # Standard Monte Carlo
        p_mean = numpy.mean(arith_payoff)
        p_std = numpy.std(arith_payoff)
        p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
        return p_mean, p_std, p_confmc

    elif control_variate == GEO_MEAN:
        dt = T / N
        sigma = V
        drift = math.exp((R - 0.5 * sigma * sigma) * dt)
        # arith_payoff = []
        # geo_payoff = []
        arith_payoff = numpy.empty(path_num)
        geo_payoff = numpy.empty(path_num)
        sigma_sqrt = sigma * math.sqrt(dt)
        exp_RT = math.exp(-R * T)
        for i in xrange(path_num):
            # spath = []
            spath = numpy.empty(N)
            former = S0
            for j in xrange(int(N)):
                # growth_factor = drift * math.exp(sigma * math.sqrt(dt) * numpy.random.normal(0, 1))
                growth_factor = drift * math.exp(sigma_sqrt * numpy.random.normal(0, 1))
                former = former * growth_factor
                # spath.append(former)
                spath[j] = former

            arith_mean = numpy.mean(spath)

            # geo_spath = map(lambda x: math.log(x), spath)
            # geo_mean = math.exp((1 / N) * sum(geo_spath))

            # geo_spath = reduce(lambda x, y: x * y, spath)
            geo_spath = numpy.prod(spath)
            geo_mean = numpy.power(geo_spath, (1 / float(len(spath))))

            if option_type == 'call':
                # arith_payoff_call = math.exp(-R * T) * max(arith_mean - K, 0)
                # geo_payoff_call = math.exp(-R * T) * max(geo_mean - geo_K, 0)
                arith_payoff_call = exp_RT * max(arith_mean - K, 0)
                geo_payoff_call = exp_RT * max(geo_mean - geo_K, 0)
                # arith_payoff.append(arith_payoff_call)
                # geo_payoff.append(geo_payoff_call)
                arith_payoff[i] = arith_payoff_call
                geo_payoff[i] = geo_payoff_call
            elif option_type == 'put':
                # arith_payoff_put = math.exp(-R * T) * max(K - arith_mean, 0)
                # geo_payoff_put = math.exp(-R * T) * max(geo_K - geo_mean, 0)
                arith_payoff_put = exp_RT * max(K - arith_mean, 0)
                geo_payoff_put = exp_RT * max(geo_K - geo_mean, 0)
                # arith_payoff.append(arith_payoff_put)
                # geo_payoff.append(geo_payoff_put)
                arith_payoff[i] = arith_payoff_put
                geo_payoff[i] = geo_payoff_put


        # Control Variate
        # covxy = numpy.mean([x * y for x, y in zip(geo_payoff, arith_payoff)]) - numpy.mean(arith_payoff) * numpy.mean(
        #     geo_payoff)
        covxy = numpy.mean(geo_payoff * arith_payoff) - numpy.mean(arith_payoff) * numpy.mean(geo_payoff)
        theta = covxy / numpy.var(geo_payoff)

        # Control Variate Version
        geo = geometric_asian_option(geo_K, T, R, V, S0, N, option_type)
        # z = [x + y for x, y in zip(arith_payoff, map(lambda x: theta * (geo - x), geo_payoff))]
        z = arith_payoff + theta * (geo - geo_payoff)
        z_mean = numpy.mean(z)
        z_std = numpy.std(z)
        z_confmc = (z_mean - 1.96 * z_std / math.sqrt(path_num), z_mean + 1.96 * z_std / math.sqrt(path_num))
        return z_mean, z_std, z_confmc

    elif control_variate == GEO_MEAN_STRIKE:
        # K = K + mean(agT)-mean(aaT)
        # mean(agT) = e^(mu*T)*ag_0 = e^(mu*T)*S0
        # mean(aaT) = 1/n * sum(S(T)) = 1/n * sum(S0*e^(rt))
        sigma = _get_geometric_sigma(V, N)
        miu = _get_geometric_miu(R, V, N, sigma)
        E_ag = S0 * math.exp(miu * T)

        dt = T / N
        E_aa = sum([math.exp(R * (i + 1) * dt) for i in xrange(int(N))]) * S0 / N
        geo_K = K + E_ag - E_aa
        return arithmetic_asian_option(K, geo_K, T, R, V, S0, N, option_type, path_num, GEO_MEAN)


def arithmetic_basket_option(S1, S2, V1, V2, R, T, K, geo_K, rou, option_type, path_num=10000,
                             control_variate='Standard'):
    ran_arg_1 = (R - 0.5 * V1 * V1) * T
    ran_arg_2 = (R - 0.5 * V2 * V2) * T
    v1_sqrt = V1 * math.sqrt(T)
    v2_sqrt = V2 * math.sqrt(T)
    exp_RT = math.exp(-R * T)
    if control_variate == STANDARD:
        # arith_basket_payoff = []
        arith_basket_payoff = numpy.empty(path_num)
        for i in xrange(path_num):
            ran1 = numpy.random.normal(0, 1)
            ran2 = rou * ran1 + math.sqrt(1 - rou * rou) * numpy.random.normal(0, 1)
            # a1 = S1 * math.exp((R - 0.5 * V1 * V1) * T + V1 * math.sqrt(T) * ran1)
            # a2 = S2 * math.exp((R - 0.5 * V2 * V2) * T + V2 * math.sqrt(T) * ran2)
            a1 = S1 * math.exp(ran_arg_1 + v1_sqrt * ran1)
            a2 = S2 * math.exp(ran_arg_2 + v2_sqrt * ran2)
            arith_basket_mean = (a1 + a2) / 2

            if option_type == 'call':
                # arith_basket_payoff_call = math.exp(-R * T) * max(arith_basket_mean - K, 0)
                arith_basket_payoff_call = exp_RT * max(arith_basket_mean - K, 0)
                arith_basket_payoff[i] = arith_basket_payoff_call
            elif option_type == 'put':
                # arith_basket_payoff_put = math.exp(-R * T) * max(K - arith_basket_mean, 0)
                arith_basket_payoff_put = exp_RT * max(K - arith_basket_mean, 0)
                arith_basket_payoff[i] = arith_basket_payoff_put


        # Standard Monte Carlo
        p_mean = numpy.mean(arith_basket_payoff)
        p_std = numpy.std(arith_basket_payoff)
        p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
        return p_mean, p_std, p_confmc

    elif control_variate == GEO_MEAN:
        # arith_basket_payoff = []
        # geo_basket_payoff = []
        arith_basket_payoff = numpy.empty(path_num)
        geo_basket_payoff = numpy.empty(path_num)
        for i in xrange(path_num):
            ran1 = numpy.random.normal(0, 1)
            ran2 = rou * ran1 + math.sqrt(1 - rou * rou) * numpy.random.normal(0, 1)
            # a1 = S1 * math.exp((R - 0.5 * V1 * V1) * T + V1 * math.sqrt(T) * ran1)
            # a2 = S2 * math.exp((R - 0.5 * V2 * V2) * T + V2 * math.sqrt(T) * ran2)
            a1 = S1 * math.exp(ran_arg_1 + v1_sqrt * ran1)
            a2 = S2 * math.exp(ran_arg_2 + v2_sqrt * ran2)
            arith_basket_mean = (a1 + a2) / 2

            # geo_basket_mean = math.exp(math.log(arith_basket_mean))
            geo_basket_mean = math.sqrt(a1 * a2)

            if option_type == 'call':
                # arith_basket_payoff_call = math.exp(-R * T) * max(arith_basket_mean - K, 0)
                arith_basket_payoff_call = exp_RT * max(arith_basket_mean - K, 0)
                arith_basket_payoff[i] = arith_basket_payoff_call
                # geo_basket_payoff_call = math.exp(-R * T) * max(geo_basket_mean - geo_K, 0)
                geo_basket_payoff_call = exp_RT * max(geo_basket_mean - geo_K, 0)
                geo_basket_payoff[i] = geo_basket_payoff_call
            elif option_type == 'put':
                # arith_basket_payoff_put = math.exp(-R * T) * max(K - arith_basket_mean, 0)
                arith_basket_payoff_put = exp_RT * max(K - arith_basket_mean, 0)
                arith_basket_payoff[i] = arith_basket_payoff_put
                # geo_basket_payoff_put = math.exp(-R * T) * max(geo_K - geo_basket_mean, 0)
                geo_basket_payoff_put = exp_RT * max(geo_K - geo_basket_mean, 0)
                geo_basket_payoff[i] = geo_basket_payoff_put

        # Control Variate
        # covxy = numpy.mean([x * y for x, y in zip(geo_basket_payoff, arith_basket_payoff)]) - numpy.mean(
        #     arith_basket_payoff) * numpy.mean(geo_basket_payoff)
        covxy = numpy.mean(geo_basket_payoff * arith_basket_payoff) - numpy.mean(arith_basket_payoff) * numpy.mean(
            geo_basket_payoff)
        theta = covxy / numpy.var(geo_basket_payoff)

        # Control Variate Version
        geo = geometric_basket_option(S1, S2, V1, V2, R, T, geo_K, rou, option_type)
        # z = [x + y for x, y in zip(arith_basket_payoff, map(lambda x: theta * (geo - x), geo_basket_payoff))]
        z = arith_basket_payoff + theta * (geo - geo_basket_payoff)
        z_mean = numpy.mean(z)
        z_std = numpy.std(z)
        z_confmc = (z_mean - 1.96 * z_std / math.sqrt(path_num), z_mean + 1.96 * z_std / math.sqrt(path_num))
        return z_mean, z_std, z_confmc

    elif control_variate == GEO_MEAN_STRIKE:
        # K = K + mean(bgT)-mean(baT)
        # mean(bgT) = e^(mu*T)*bg_0 = e^(mu*T)*S0
        # mean(baT) = 1/n * sum(S(T)) = 1/n * sum(S0*e^(rt))
        bg_0 = math.sqrt(S1 * S2)
        sigma_bg = math.sqrt(V1 * V1 + V1 * V2 * rou + V2 * V1 * rou + V2 * V2) / 2
        miu_bg = R - 0.5 * (V1 * V1 + V2 * V2) / 2 + 0.5 * sigma_bg * sigma_bg

        E_bg = bg_0 * math.exp(miu_bg * T)
        E_ba = (S1 * math.exp(R * T) + S2 * math.exp(R * T)) / 2
        geo_K = K + E_bg - E_ba
        return arithmetic_basket_option(S1, S2, V1, V2, R, T, K, geo_K, rou, option_type, path_num, GEO_MEAN)


if __name__ == '__main__':
    #print"S=100,K=100,t=0,T=0.5,v=20%,and r=1%."
	pass
    # S = S0 = S1 = S2 = 100.0
#     T = 3.0
#     R = 0.05
#     V = V1 = V2 = 0.3
#     K = 100.0
#     n = 50.0
#     rou = 0.5
#     m = 10000
# 
#     call_option_price = bs(100, 100, 0.5, 0.2, 0.01, 'call')
#     put_option_price = bs(100, 100, 0.5, 0.2, 0.01, 'put')
# 
#     call_option_asian = geometric_asian_option(K, T, R, V, S0, n, 'call')
#     put_option_asian = geometric_asian_option(K, T, R, V, S0, n, 'put')
# 
#     call_option_basket = geometric_basket_option(S1, S2, V1, V2, R, T, K, rou, 'call')
#     put_option_basket = geometric_basket_option(S1, S2, V1, V2, R, T, K, rou, 'put')
# 
#     print "European call/put options:"
#     print "call:", call_option_price
#     print "put", put_option_price
# 
#     print "geometric Asian call/put options:"
#     print "1call:", call_option_asian
#     print "1put", put_option_asian
# 
#     print "geometric_basket_asian_option:"
#     print "2call:", call_option_basket
#     print "2put", put_option_basket
# 
#     MC_call_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'call', m, 'Standard')
#     MC_put_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'put', m, 'Standard')
#     print "arithmetic_asian_option-Standard Monte Carlo:"
#     print "1call:", MC_call_arithmetic_asian
#     print "1put", MC_put_arithmetic_asian
# 
#     MC_call_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'call', m, 'Standard')
#     MC_put_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'put', m, 'Standard')
#     print "arithmetic_basket_asian_option-Standard Monte Carlo"
#     print "2call:", MC_call_arithmetic_basket
#     print "2put", MC_put_arithmetic_basket
# 
#     MCG_call_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'call', m, 'Geometric mean Asian')
#     MCG_put_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'put', m, 'Geometric mean Asian')
#     print "arithmetic_asian_option-Control Variate:"
#     print "1call:", MCG_call_arithmetic_asian
#     print "1put", MCG_put_arithmetic_asian
# 
#     MCG_call_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'call', m,
#                                                           'Geometric mean Asian')
#     MCG_put_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'put', m,
#                                                          'Geometric mean Asian')
#     print "arithmetic_basket_asian_option-Control Variate"
#     print "2call:", MCG_call_arithmetic_basket
#     print "2put", MCG_put_arithmetic_basket
# 
#     astrike_call_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'call', m,
#                                                             'Geometric mean Asian with adjusted strike')
#     astrike_put_arithmetic_asian = arithmetic_asian_option(K, K, T, R, V, S0, n, 'put', m,
#                                                            'Geometric mean Asian with adjusted strike')
#     print "arithmetic_asian_option-Control Variate:"
#     print "1call:", astrike_call_arithmetic_asian
#     print "1put", astrike_put_arithmetic_asian
# 
#     astrike_call_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'call', m,
#                                                               'Geometric mean Asian with adjusted strike')
#     astrike_put_arithmetic_basket = arithmetic_basket_option(S1, S2, V1, V2, R, T, K, K, rou, 'put', m,
#                                                              'Geometric mean Asian with adjusted strike')
#     print "arithmetic_basket_asian_option-Control Variate"
#     print "2call:", astrike_call_arithmetic_basket
#     print "2put", astrike_put_arithmetic_basket
# 
#     # print "1call", call_option_asian, ":::::::::", MC_call_arithmetic_asian, "::::::::::", MCG_call_arithmetic_asian, ":::::::::::", astrike_call_arithmetic_asian
#     print "................."
#     print MC_call_arithmetic_basket, MCG_call_arithmetic_basket, astrike_call_arithmetic_basket
#     print MC_put_arithmetic_basket, MCG_put_arithmetic_basket, astrike_put_arithmetic_basket
# 
#     print "--------------------"
# 
#     print "arithmetic_asian_option-Standard Monte Carlo:"
#     print "1call:", MC_call_arithmetic_asian
#     print "1put", MC_put_arithmetic_asian
# 
#     print "arithmetic_asian_option-Control Variate:"
#     print "1call:", MCG_call_arithmetic_asian
#     print "1put", MCG_put_arithmetic_asian
# 
#     print "arithmetic_asian_option-Control Variate:"
#     print "1call:", astrike_call_arithmetic_asian
#     print "1put", astrike_put_arithmetic_asian

