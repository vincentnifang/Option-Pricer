__author__ = 'vincent'
from scipy import stats
import pyopencl as cl
import numpy, math, time
import Quasi_Monte_Carlo as quasi

STANDARD = 'Standard'
GEO_MEAN = 'Geometric mean Asian'
GEO_MEAN_STRIKE = 'Geometric mean Asian with adjusted strike'


def geometric_basket_option(S1, S2, V1, V2, R, T, K, rou, option_type):
    bg_0 = math.sqrt(S1 * S2)
    sigma_bg = math.sqrt(V1 * V1 + V1 * V2 * rou + V2 * V1 * rou + V2 * V2) / 2
    miu_bg = R - 0.5 * (V1 * V1 + V2 * V2) / 2 + 0.5 * sigma_bg * sigma_bg

    d1 = (math.log(bg_0 / K) + (miu_bg + 0.5 * sigma_bg * sigma_bg) * T) / (sigma_bg * math.sqrt(T))
    d2 = d1 - sigma_bg * math.sqrt(T)

    if option_type == 1.0:
        return math.exp(-R * T) * (bg_0 * math.exp(miu_bg * T) * stats.norm.cdf(d1) - K * stats.norm.cdf(d2))
    elif option_type == 2.0:
        return math.exp(-R * T) * (K * stats.norm.cdf(-d2) - bg_0 * math.exp(miu_bg * T) * stats.norm.cdf(-d1))
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
    if option_type == 1.0:
        return math.exp(-R * T) * (S0 * math.exp(miu * T) * stats.norm.cdf(d1) - K * stats.norm.cdf(d2))
    elif option_type == 2.0:
        return math.exp(-R * T) * (K * stats.norm.cdf(-d2) - S0 * math.exp(miu * T) * stats.norm.cdf(-d1))
    return None


def GPU_arithmetic_basket_option(S1, S2, V1, V2, R, T, K, geo_K, rou, option_type, path_num=10000,
                                 control_variate='Standard', Quasi=True):
    if control_variate == STANDARD:
        cntxt = cl.create_some_context()
        #now create a command queue in the context
        queue = cl.CommandQueue(cntxt)
        # create some data array to give as input to Kernel and get output
        # rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)
        # rand2 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)
        if Quasi == True:
            rand1 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num), 2.0), dtype=numpy.float32)
            rand2 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num), 2.0), dtype=numpy.float32)
        else:
            rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)
            rand2 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)

        arith_basket_payoff = numpy.empty(rand1.shape, dtype=numpy.float32)
        # create the buffers to hold the values of the input
        rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
        rand2_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand2)
        # create output buffer
        arith_basket_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, arith_basket_payoff.nbytes)


        # Kernel Program
        code = """
        __kernel void standard_arithmetic_basket_option(__global float* num1, __global float* num2, __global float* arith_basket_payoff ,float S1, float S2, float V1, float V2, float R, float K, float T, float rou, float option_type)
        {
            int i = get_global_id(0);
            float ran1 = num1[i];
            float ran2 = rou * ran1 + sqrt(1 - rou * rou) * num2[i];
            float a1 = S1 * exp((R - 0.5 * V1 * V1) * T + V1 * sqrt(T) * ran1);
            float a2 = S2 * exp((R - 0.5 * V2 * V2) * T + V2 * sqrt(T) * ran2);
            float arith_basket_mean = (a1 + a2) / 2;

            if isequal(option_type, 1.0){
                float arith_basket_payoff_call = exp(-R * T) * fmax(arith_basket_mean - K, 0);
                arith_basket_payoff[i] = arith_basket_payoff_call;
            }
            else if isequal(option_type, 2.0){
                float arith_basket_payoff_put = exp(-R * T) * fmax(K - arith_basket_mean, 0);
                arith_basket_payoff[i] = arith_basket_payoff_put;
            }

        }
        """
        # build the Kernel
        bld = cl.Program(cntxt, code).build()
        S1 = numpy.float32(S1)
        S2 = numpy.float32(S2)
        V1 = numpy.float32(V1)
        V2 = numpy.float32(V2)
        R = numpy.float32(R)
        K = numpy.float32(K)
        T = numpy.float32(T)
        rou = numpy.float32(rou)
        option_type = numpy.float32(option_type)

        kernelargs = (S1, S2, V1, V2, R, K, T, rou, option_type)
        # Kernel is now launched
        launch = bld.standard_arithmetic_basket_option(queue, rand1.shape, None, rand1_buf, rand2_buf,
                                                       arith_basket_payoff_buf, *(kernelargs))
        # wait till the process completes
        launch.wait()
        cl.enqueue_read_buffer(queue, arith_basket_payoff_buf, arith_basket_payoff).wait()
        # print the output
        p_mean = numpy.mean(arith_basket_payoff)
        p_std = numpy.std(arith_basket_payoff)
        p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
        return p_mean, p_std, p_confmc

    elif control_variate == GEO_MEAN:
        cntxt = cl.create_some_context()
        #now create a command queue in the context
        queue = cl.CommandQueue(cntxt)
        # create some data array to give as input to Kernel and get output
        # rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)
        # rand2 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)

        if Quasi == True:
            rand1 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num), 2.0), dtype=numpy.float32)
            rand2 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num), 2.0), dtype=numpy.float32)
        else:
            rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)
            rand2 = numpy.array(numpy.random.normal(0, 1, (path_num, 1)), dtype=numpy.float32)

        arith_basket_payoff = numpy.empty(rand1.shape, dtype=numpy.float32)
        geo_basket_payoff = numpy.empty(rand1.shape, dtype=numpy.float32)
        # create the buffers to hold the values of the input
        rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
        rand2_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand2)
        # create output buffer
        arith_basket_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, arith_basket_payoff.nbytes)
        geo_basket_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, geo_basket_payoff.nbytes)

        # Kernel Program
        code = """
        __kernel void geo_mean_arithmetic_basket_option(__global float* num1, __global float* num2, __global float* out,__global float* geo, float S1, float S2, float V1, float V2, float R, float K,float geo_K, float T, float rou, float option_type)
        {
            int i = get_global_id(0);
            float ran1 = num1[i];
            float ran2 = rou * ran1 + sqrt(1 - rou * rou) * num2[i];
            float a1 = S1 * exp((R - 0.5 * V1 * V1) * T + V1 * sqrt(T) * ran1);
            float a2 = S2 * exp((R - 0.5 * V2 * V2) * T + V2 * sqrt(T) * ran2);
            float arith_basket_mean = (a1 + a2) / 2;

            float geo_basket_mean = sqrt(a1 * a2);

            if isequal(option_type, 1.0){
                float arith_basket_payoff_call = exp(-R * T) * fmax(arith_basket_mean - K, 0);
                out[i] = arith_basket_payoff_call;
                float geo_basket_payoff_call = exp(-R * T) * fmax(geo_basket_mean - geo_K, 0);
                geo[i] = geo_basket_payoff_call;
            }
            else if isequal(option_type, 2.0){
                float arith_basket_payoff_put = exp(-R * T) * fmax(K - arith_basket_mean, 0);
                out[i] = arith_basket_payoff_put;
                float geo_basket_payoff_put = exp(-R * T) * fmax(geo_K - geo_basket_mean, 0);
                geo[i] = geo_basket_payoff_put;
            }

        }
        """
        # build the Kernel
        bld = cl.Program(cntxt, code).build()
        S1 = numpy.float32(S1)
        S2 = numpy.float32(S2)
        V1 = numpy.float32(V1)
        V2 = numpy.float32(V2)
        R = numpy.float32(R)
        K = numpy.float32(K)
        T = numpy.float32(T)
        rou = numpy.float32(rou)
        option_type = numpy.float32(option_type)
        geo_K = numpy.float32(geo_K)

        kernelargs = (S1, S2, V1, V2, R, K, geo_K, T, rou, option_type)
        # Kernel is now launched
        launch = bld.geo_mean_arithmetic_basket_option(queue, rand1.shape, None, rand1_buf, rand2_buf,
                                                       arith_basket_payoff_buf, geo_basket_payoff_buf, *(kernelargs))
        # wait till the process completes
        launch.wait()
        cl.enqueue_read_buffer(queue, arith_basket_payoff_buf, arith_basket_payoff).wait()
        cl.enqueue_read_buffer(queue, geo_basket_payoff_buf, geo_basket_payoff).wait()
        # print the output

        # Control Variate
        covxy = numpy.mean(geo_basket_payoff * arith_basket_payoff) - numpy.mean(
            arith_basket_payoff) * numpy.mean(geo_basket_payoff)
        theta = covxy / numpy.var(geo_basket_payoff)

        # Control Variate Version
        geo = geometric_basket_option(S1, S2, V1, V2, R, T, geo_K, rou, option_type)
        z = arith_basket_payoff + theta * (geo - geo_basket_payoff)
        # z = [x + y for x, y in zip(arith_basket_payoff, map(lambda x: theta * (geo - x), geo_basket_payoff))]
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
        return GPU_arithmetic_basket_option(S1, S2, V1, V2, R, T, K, geo_K, rou, option_type, path_num, GEO_MEAN)


def GPU_arithmetic_asian_option(K, geo_K, T, R, V, S0, N, option_type, path_num=10000, control_variate='Standard',
                                Quasi=True):
    if control_variate == STANDARD:
        dt = T / N
        sigma = V
        drift = math.exp((R - 0.5 * sigma * sigma) * dt)
        sigma_sqrt = sigma * math.sqrt(dt)
        exp_RT = math.exp(-R * T)

        cntxt = cl.create_some_context()
        #now create a command queue in the context
        queue = cl.CommandQueue(cntxt)
        # create some data array to give as input to Kernel and get output
        # rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)
        if Quasi == True:
            # rand1 = numpy.array(quasi.quasi_normal_random(int(path_num * N), 2.0), dtype=numpy.float32)
            rand1 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num * N), 2.0), dtype=numpy.float32)
        else:
            rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)

        arith_asian_payoff = numpy.empty((path_num, 1), dtype=numpy.float32)
        # create the buffers to hold the values of the input
        rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
        # create output buffer
        arith_asian_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, arith_asian_payoff.nbytes)

        # Kernel Program
        code = """
        __kernel void standard_arithmetic_asian_option(__global float* num1, __global float* arith_asian_payoff,float N, float K, float S0, float sigma_sqrt, float drift, float exp_RT,float option_type)
        {
            int i = get_global_id(0);
            float former = S0;
            float arith_asian_all = 0.0;
            for(int j=0;j<(int)N;j++){
                float rand = num1[i*(int)N+j];
                float growth_factor = drift * exp(sigma_sqrt * rand);
                former = former * growth_factor;
                arith_asian_all = arith_asian_all + former;
            }
            float arith_asian_mean = arith_asian_all / N;

            if isequal(option_type, 1.0){
                arith_asian_payoff[i] = exp_RT * fmax(arith_asian_mean - K, 0);
            }
            else if isequal(option_type, 2.0){
                arith_asian_payoff[i] = exp_RT * fmax(K - arith_asian_mean, 0);
            }

        }
        """
        # build the Kernel
        bld = cl.Program(cntxt, code).build()
        N = numpy.float32(N)
        K = numpy.float32(K)
        S0 = numpy.float32(S0)
        sigma_sqrt = numpy.float32(sigma_sqrt)
        drift = numpy.float32(drift)
        exp_RT = numpy.float32(exp_RT)
        option_type = numpy.float32(option_type)

        kernelargs = (N, K, S0, sigma_sqrt, drift, exp_RT, option_type)
        # Kernel is now launched
        # print rand1.shape
        launch = bld.standard_arithmetic_asian_option(queue, (path_num, 1), None, rand1_buf,
                                                      arith_asian_payoff_buf, *(kernelargs))
        # wait till the process completes
        launch.wait()
        cl.enqueue_read_buffer(queue, arith_asian_payoff_buf, arith_asian_payoff).wait()

        # print the output
        p_mean = numpy.mean(arith_asian_payoff)
        p_std = numpy.std(arith_asian_payoff)
        p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
        return p_mean, p_std, p_confmc
    elif control_variate == GEO_MEAN:
        dt = T / N
        sigma = V
        drift = math.exp((R - 0.5 * sigma * sigma) * dt)
        sigma_sqrt = sigma * math.sqrt(dt)
        exp_RT = math.exp(-R * T)

        cntxt = cl.create_some_context()
        #now create a command queue in the context
        queue = cl.CommandQueue(cntxt)
        # create some data array to give as input to Kernel and get output
        # rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)
        if Quasi == True:
            rand1 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num * N), 2.0), dtype=numpy.float32)
        else:
            rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)

        arith_payoff = numpy.empty((path_num, 1), dtype=numpy.float32)
        geo_payoff = numpy.empty((path_num, 1), dtype=numpy.float32)
        # create the buffers to hold the values of the input
        rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                     cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
        # create output buffer
        arith_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, arith_payoff.nbytes)
        geo_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, geo_payoff.nbytes)


        # Kernel Program
        code = """
        __kernel void geo_mean_arithmetic_asian_option(__global float* num1, __global float* arith_payoff,__global float* geo_payoff, float N, float K,float geo_K, float S0, float sigma_sqrt, float drift, float exp_RT,float option_type)
        {
            int i = get_global_id(0);
            float former = S0;
            float arith_asian_all = 0.0;
            float geo_asian_all = 0.0;
            for(int j=0;j<(int)N;j++){
                float rand = num1[i*(int)N+j];
                float growth_factor = drift * exp(sigma_sqrt * rand);
                former = former * growth_factor;
                arith_asian_all = arith_asian_all + former;
                geo_asian_all = geo_asian_all + log(former);
            }
            float arith_asian_mean = arith_asian_all / N;
            float geo_asian_mean = exp(geo_asian_all / N);

            if isequal(option_type, 1.0){
                arith_payoff[i] = exp_RT * fmax(arith_asian_mean - K, 0);
                geo_payoff[i] = exp_RT * fmax(geo_asian_mean - geo_K, 0);
            }
            else if isequal(option_type, 2.0){
                arith_payoff[i] = exp_RT * fmax(K - arith_asian_mean, 0);
                geo_payoff[i] = exp_RT * fmax(geo_K - geo_asian_mean, 0);
            }

        }
        """
        # build the Kernel
        bld = cl.Program(cntxt, code).build()
        N = numpy.float32(N)
        K = numpy.float32(K)
        geo_K = numpy.float32(geo_K)
        S0 = numpy.float32(S0)
        sigma_sqrt = numpy.float32(sigma_sqrt)
        drift = numpy.float32(drift)
        exp_RT = numpy.float32(exp_RT)
        option_type = numpy.float32(option_type)

        kernelargs = (N, K, geo_K, S0, sigma_sqrt, drift, exp_RT, option_type)
        # Kernel is now launched
        # print rand1.shape
        launch = bld.geo_mean_arithmetic_asian_option(queue, (path_num, 1), None, rand1_buf,
                                                      arith_payoff_buf, geo_payoff_buf, *(kernelargs))
        # wait till the process completes
        launch.wait()
        cl.enqueue_read_buffer(queue, arith_payoff_buf, arith_payoff).wait()
        cl.enqueue_read_buffer(queue, geo_payoff_buf, geo_payoff).wait()
        # print the output
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
        return GPU_arithmetic_asian_option(K, geo_K, T, R, V, S0, N, option_type, path_num, GEO_MEAN)


def european_option(K, T, R, V, S0, N, option_type, path_num=10000):
    dt = T / N
    sigma = V
    drift = math.exp((R - 0.5 * sigma * sigma) * dt)
    sigma_sqrt = sigma * math.sqrt(dt)
    exp_RT = math.exp(-R * T)
    european_payoff = []
    for i in xrange(path_num):
        former = S0
        for j in xrange(int(N)):
            former = former * drift * math.exp(sigma_sqrt * numpy.random.normal(0, 1))
        european_option = former

        if option_type == 'call':
            european_payoff_call = exp_RT * max(european_option - K, 0)
            european_payoff.append(european_payoff_call)
        elif option_type == 'put':
            european_payoff_put = exp_RT * max(K - european_option, 0)
            european_payoff.append(european_payoff_put)

    # Standard Monte Carlo
    p_mean = numpy.mean(european_payoff)
    p_std = numpy.std(european_payoff)
    p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
    return p_mean, p_std, p_confmc

def GPU_european_option(K, T, R, V, S0, N, option_type, path_num=10000, Quasi=True):
    dt = T / N
    sigma = V
    drift = math.exp((R - 0.5 * sigma * sigma) * dt)
    sigma_sqrt = sigma * math.sqrt(dt)
    exp_RT = math.exp(-R * T)
    cntxt = cl.create_some_context()
    #now create a command queue in the context
    queue = cl.CommandQueue(cntxt)
    # create some data array to give as input to Kernel and get output
    # rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)
    if Quasi == True:
        # rand1 = numpy.array(quasi.quasi_normal_random(int(path_num * N), 2.0), dtype=numpy.float32)
        rand1 = numpy.array(quasi.GPU_quasi_normal_random(int(path_num * N), 2.0), dtype=numpy.float32)
    else:
        rand1 = numpy.array(numpy.random.normal(0, 1, (path_num, N)), dtype=numpy.float32)

    european_payoff = numpy.empty((path_num, 1), dtype=numpy.float32)
    # create the buffers to hold the values of the input
    rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                 cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
    # create output buffer
    european_payoff_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, european_payoff.nbytes)

    # Kernel Program
    code = """
    __kernel void european_option(__global float* num1, __global float* european_payoff,float N, float K, float S0, float sigma_sqrt, float drift, float exp_RT,float option_type)
    {
        int i = get_global_id(0);
        float former = S0;
        for(int j=0;j<(int)N;j++){
            float rand = num1[i*(int)N+j];
            former = former * drift * exp(sigma_sqrt * rand);
        }
        float european = former;
        if isequal(option_type, 1.0){
            european_payoff[i] = exp_RT * fmax(european - K, 0);
        }
        else if isequal(option_type, 2.0){
            european_payoff[i] = exp_RT * fmax(K - european, 0);
        }

    }
    """
    # build the Kernel
    bld = cl.Program(cntxt, code).build()
    N = numpy.float32(N)
    K = numpy.float32(K)
    S0 = numpy.float32(S0)
    sigma_sqrt = numpy.float32(sigma_sqrt)
    drift = numpy.float32(drift)
    exp_RT = numpy.float32(exp_RT)
    option_type = numpy.float32(option_type)

    kernelargs = (N, K, S0, sigma_sqrt, drift, exp_RT, option_type)
    # Kernel is now launched
    # print rand1.shape
    launch = bld.european_option(queue, (path_num, 1), None, rand1_buf,
                                                  european_payoff_buf, *(kernelargs))
    # wait till the process completes
    launch.wait()
    cl.enqueue_read_buffer(queue, european_payoff_buf, european_payoff).wait()

    # Standard Monte Carlo
    p_mean = numpy.mean(european_payoff)
    p_std = numpy.std(european_payoff)
    p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
    return p_mean, p_std, p_confmc


if __name__ == '__main__':
    #print"S=100,K=100,t=0,T=0.5,v=20%,and r=1%."
    import time

    s = time.time()

    S = S0 = S1 = S2 = 100.0
    T = 3.0
    R = 0.05
    V = V1 = V2 = 0.3
    K = 100.0
    n = 50.0
    rou = 0.5
    m = 10000

    print GPU_arithmetic_asian_option(K, K, T, R, V, S0, n, 1.0, path_num=100, control_variate=STANDARD, Quasi=True)

    import project

    print project.arithmetic_asian_option(K, K, T, R, V, S0, n, 'call', path_num=100, control_variate=STANDARD)

    e = time.time()
    print "use", e - s

    s = time.time()
    print european_option(K, T, R, V, S0, n, 'call', path_num=100000)
    e = time.time()
    print "use", e - s

    print project.bs(S0, K, T, V, R, 'call')
    s = time.time()
    print GPU_european_option(K, T, R, V, S0, n, 1.0, path_num=100000, Quasi=True)
    e = time.time()
    print "use", e - s