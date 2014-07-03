__author__ = 'vincent'
import math,numpy,random,time

STANDARD = 'Standard'
GEO_MEAN = 'Geometric mean Asian'
GEO_MEAN_STRIKE = 'Geometric mean Asian with adjusted strike'


a0 = 2.50662823884
b0 = -8.47351093090
a1 = -18.61500062529
b1 = 23.08336743743
a2 = 41.39119773534
b2 = -21.06224101826
a3 = -25.44106049637
b3 = 3.13082909833
c0 = 0.3374754822726147
c5 = 0.0003951896511919
c1 = 0.9761690190917186
c6 = 0.0000321767881768
c2 = 0.1607979714918209
c7 = 0.0000002888167364
c3 = 0.0276438810333863
c8 = 0.0000003960315187
c4 = 0.0038405729373609


# Beasley_Springer_Moro

def bs_moro(u):
    y = u - 0.5
    if abs(y) < 0.42:
        r = y * y
        x = y *  (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r +1)
    else:
        r = u
        if y > 0:
            r = 1 - u
        r = math.log(-math.log(r))
        x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4+ r * (c5 + r * (c6 + r * (c7 + r * c8)))))))
        if y < 0:
            x = -x
    return x

# Halton sequence

def halton(index, base):
    result = 0.0
    f = 1 / base
    i = index
    while i > 0:
        result = result + f * (i % base)
        i = numpy.floor(i / base)
        f = f / base
    return result

#
# FUNCTION (index, base)
#    BEGIN
#        result = 0;
#        f = 1 / base;
#        i = index;
#        WHILE (i > 0)
#        BEGIN
#            result = result + f * (i % base);
#            i = FLOOR(i / base);
#            f = f / base;
#        END
#        RETURN result;
#    END

def quasi_normal_random(N,base=2.0):
    result = []
    for i in xrange(N):
        result.append(bs_moro(halton(i+1, base)))
    random.shuffle(result)
    return result

def GPU_quasi_normal_random(N,base=2.0):
    import pyopencl as cl
    platform = cl.get_platforms()
    my_gpu_devices = platform[0].get_devices(device_type=cl.device_type.GPU)
    cntxt = cl.Context(devices=my_gpu_devices)

    # cntxt = cl.create_some_context()
    #now create a command queue in the context
    queue = cl.CommandQueue(cntxt)
    # create some data array to give as input to Kernel and get output
    rand1 = numpy.empty((N, 1), dtype=numpy.float32)
    out = numpy.empty((N, 1), dtype=numpy.float32)
    # create the buffers to hold the values of the input
    rand1_buf = cl.Buffer(cntxt, cl.mem_flags.READ_ONLY |
                                 cl.mem_flags.COPY_HOST_PTR, hostbuf=rand1)
    # create output buffer
    out_buf = cl.Buffer(cntxt, cl.mem_flags.WRITE_ONLY, out.nbytes)

    # Kernel Program
    code = """
    __kernel void GPU_quasi(__global float* num1, __global float* out,float base)
    {
        int i = get_global_id(0);
        float result = 0.0;
        float f = 1 / base;
        float index = i + 1.0;
        while(index > 0) {
            result = result + f * ((int)index % (int)base);
            index = floor(index / base);
            f = f / base;
        }


        float a0 = 2.50662823884;
        float b0 = -8.47351093090;
        float a1 = -18.61500062529;
        float b1 = 23.08336743743;
        float a2 = 41.39119773534;
        float b2 = -21.06224101826;
        float a3 = -25.44106049637;
        float b3 = 3.13082909833;
        float c0 = 0.3374754822726147;
        float c5 = 0.0003951896511919;
        float c1 = 0.9761690190917186;
        float c6 = 0.0000321767881768;
        float c2 = 0.1607979714918209;
        float c7 = 0.0000002888167364;
        float c3 = 0.0276438810333863;
        float c8 = 0.0000003960315187;
        float c4 = 0.0038405729373609;

        float y = result - 0.5;
        float r = 0.0;
        float x = 0.0;
        if (fabs(y) < 0.42) {
            r = y * y;
            x = y *  (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r +1);
        }
        else{
            r = result;
            if (y > 0) {
                r = 1 - result;
            }
            r = log(-log(r));
            x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4+ r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
            if (y < 0){
                x = -x;
            }
        }
        out[i] = x;

    }
    """
    # build the Kernel
    bld = cl.Program(cntxt, code).build()
    base = numpy.float32(base)


    # Kernel is now launched
    # print rand1.shape
    launch = bld.GPU_quasi(queue, (N, 1), None, rand1_buf,out_buf, base)
    # wait till the process completes
    launch.wait()
    cl.enqueue_read_buffer(queue, out_buf, out).wait()
    random.shuffle(out)
    return out

if __name__ == "__main__":
    N = 10000
    base = 2.0
    s = time.time()
    print numpy.mean(GPU_quasi_normal_random(N,base))
    e = time.time()
    print e-s

    s = time.time()
    print numpy.mean(quasi_normal_random(N,base))
    e = time.time()
    print e-s


    # S = S0 = S1 = S2 = 100.0
    # T = 3.0
    # R = 0.05
    # V = V1 = V2 = 0.3
    # K = 100.0
    # n = 50.0
    # rou = 0.5
    # m = 10000
    # geo_K = K
    # option_type = 'call'
    # control_variate = GEO_MEAN
    # path_num = 10000
    #
    #
    # print "using Quasi MC","path number=",path_num
    # start = time.time()
    #
    # random_list1 = quasi_normal_random(path_num,2.0)
    # random_list2 = quasi_normal_random(path_num,2.0)
    #
    #
    # # random_list1 = GPU_QMC.generate_GPU_QMC_random(path_num)
    # # random_list2 = GPU_QMC.generate_GPU_QMC_random(path_num)
    #
    # end = time.time()
    # print "generating random number:",end-start
    #
    # arith_basket_payoff = []
    # for i in xrange(path_num):
    #     ran1 = random_list1[i]
    #     ran2 = rou * ran1 + math.sqrt(1 - rou * rou) * random_list2[i]
    #     a1 = S1 * math.exp((R - 0.5 * V1 * V1) * T + V1 * math.sqrt(T) * ran1)
    #     a2 = S2 * math.exp((R - 0.5 * V2 * V2) * T + V2 * math.sqrt(T) * ran2)
    #     arith_basket_mean = (a1 + a2) / 2
    #
    #     if option_type == 'call':
    #         arith_basket_payoff_call = math.exp(-R * T) * max(arith_basket_mean - K, 0)
    #         arith_basket_payoff.append(arith_basket_payoff_call)
    #     elif option_type == 'put':
    #         arith_basket_payoff_put = math.exp(-R * T) * max(K - arith_basket_mean, 0)
    #         arith_basket_payoff.append(arith_basket_payoff_put)
    #
    #
    # # Standard Monte Carlo
    # p_mean = numpy.mean(arith_basket_payoff)
    # p_std = numpy.std(arith_basket_payoff)
    # p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
    #
    #
    # end = time.time()
    # print "time", end - start
    # print p_mean, p_std, p_confmc
    #
    # print "--------------------------------"
    # print "using standard MC","path number=",path_num
    #
    # start = time.time()
    #
    #
    # arith_basket_payoff = []
    # for i in xrange(path_num):
    #     ran1 = numpy.random.normal(0, 1)
    #     ran2 = rou * ran1 + math.sqrt(1 - rou * rou) * numpy.random.normal(0, 1)
    #     a1 = S1 * math.exp((R - 0.5 * V1 * V1) * T + V1 * math.sqrt(T) * ran1)
    #     a2 = S2 * math.exp((R - 0.5 * V2 * V2) * T + V2 * math.sqrt(T) * ran2)
    #     arith_basket_mean = (a1 + a2) / 2
    #
    #     if option_type == 'call':
    #         arith_basket_payoff_call = math.exp(-R * T) * max(arith_basket_mean - K, 0)
    #         arith_basket_payoff.append(arith_basket_payoff_call)
    #     elif option_type == 'put':
    #         arith_basket_payoff_put = math.exp(-R * T) * max(K - arith_basket_mean, 0)
    #         arith_basket_payoff.append(arith_basket_payoff_put)
    #
    #
    # # Standard Monte Carlo
    # p_mean = numpy.mean(arith_basket_payoff)
    # p_std = numpy.std(arith_basket_payoff)
    # p_confmc = (p_mean - 1.96 * p_std / math.sqrt(path_num), p_mean + 1.96 * p_std / math.sqrt(path_num))
    #
    # end = time.time()
    # print "time", end - start
    # print p_mean, p_std, p_confmc
    #
    # print "--------------------------------"