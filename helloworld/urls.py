from django.conf.urls import patterns, include, url

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'helloworld.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

   # url(r'^admin/', include(admin.site.urls)),
	('^hello/$', 'hello.views.hello'),
    ('^time/$', 'hello.views.current_datetime'),
    ('^option_pricer/$', 'hello.views.option_pricer'),
    ('^bs_euro/$','hello.views.bs_euro'),
    ('^get_geometric_asian_option/$','hello.views.get_geometric_asian_option'),
    ('^get_arithmetic_asian_option/$','hello.views.get_arithmetic_asian_option'),
    ('^get_geometric_basket_option/$','hello.views.get_geometric_basket_option'),
    ('^get_arithmetic_basket_option/$','hello.views.get_arithmetic_basket_option'),
    ('^option_pricer_premium/$', 'hello.views.goto_premium'),
    ('^get_standardMC_european_option/$','hello.views.get_standardMC_european_option'),
    ('^get_GPU_european_option/$','hello.views.get_GPU_european_option'),
    ('^get_GPU_arithmetic_basket_option/$','hello.views.get_GPU_arithmetic_basket_option'),
    ('^get_GPU_arithmetic_asian_option/$','hello.views.get_GPU_arithmetic_asian_option'),
)
