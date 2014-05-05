from django.shortcuts import render_to_response

import project
import pdb


# Create your views here.
from django.http import HttpResponse
import datetime


def hello(request):
    return HttpResponse("Hello world")


def current_datetime(request):
    now = datetime.datetime.now()
    html = "<html><body>It is now %s.</body></html>" % now
    return HttpResponse(html)


#
# def search(request):
#     if 'q' in request.GET:
#         message = 'You searched for: %r' % request.GET['q']
#     else:
#         message = 'You submitted an empty form.'
#     return HttpResponse(message)


def option_pricer(request):
    return render_to_response('index.html')


def bs_euro(request):
    spot_price = request.GET.get("spot_price")
    volatility = request.GET.get("volatility")
    rate = request.GET.get("rate")
    maturity = request.GET.get("maturity")
    strike_price = request.GET.get("strike_price")
    option_type = request.GET.get("option_type")
    price = project.bs(spot_price, strike_price, maturity, volatility, rate, option_type)
    return HttpResponse(price)


def get_geometric_asian_option(request):
    # pdb.set_trace()
    spot_price = float(request.GET.get("spot_price"))
    volatility = float(request.GET.get("volatility"))
    rate = float(request.GET.get("rate"))
    maturity = float(request.GET.get("maturity"))
    strike_price = float(request.GET.get("strike_price"))
    option_type = request.GET.get("option_type")
    observation_num = float(request.GET.get("observation_num"))
    #     price=project.geometric_asian_option(K, T, R, V, S0, N, option_type)

    price = project.geometric_asian_option(strike_price, maturity, rate, volatility, spot_price, observation_num,
                                           option_type)
    return HttpResponse(price)



def get_geometric_basket_option(request):
    S1 = float(request.GET.get("S1"))
    S2 = float(request.GET.get("S2"))
    V1 = float(request.GET.get("V1"))
    V2 = float(request.GET.get("V2"))
    R = float(request.GET.get("R"))
    T = float(request.GET.get("T"))
    K = float(request.GET.get("K"))
    rou = float(request.GET.get("rou"))
    option_type = request.GET.get("option_type")
    price = project.geometric_basket_option(S1, S2, V1, V2, R, T, K, rou, option_type)
    return HttpResponse(price)

def get_arithmetic_asian_option(request):
    S0 = float(request.GET.get("spot_price"))
    V = float(request.GET.get("volatility"))
    R = float(request.GET.get("rate"))
    T = float(request.GET.get("maturity"))
    K = float(request.GET.get("strike_price"))
    geo_K = K
    option_type = request.GET.get("option_type")
    N = float(request.GET.get("observation_num"))
    path_num = int(request.GET.get("path_num"))
    control_variate = request.GET.get("control_variate")
    price, std, confcv = project.arithmetic_asian_option(K, geo_K, T, R, V, S0, N, option_type, path_num,
                                                         control_variate)
    return HttpResponse(str(price)+':'+str(confcv))



def get_arithmetic_basket_option(request):
    S1 = float(request.GET.get("S1"))
    S2 = float(request.GET.get("S2"))
    V1 = float(request.GET.get("V1"))
    V2 = float(request.GET.get("V2"))
    R = float(request.GET.get("R"))
    T = float(request.GET.get("T"))
    K = float(request.GET.get("K"))
    geo_K = K
    rou = float(request.GET.get("rou"))
    option_type = request.GET.get("option_type")
    path_num = int(request.GET.get("path_num"))
    control_variate = request.GET.get("control_variate")
    price, std, confcv = project.arithmetic_basket_option(S1, S2, V1, V2, R, T, K, geo_K, rou, option_type, path_num, control_variate)
    return HttpResponse(str(price)+':'+str(confcv))

