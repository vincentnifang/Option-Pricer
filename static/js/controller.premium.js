var MC_euro_option_url = "/get_standardMC_european_option";
var GPU_euro_option_url = "/get_GPU_european_option";
var GPU_arith_option_url = "/get_GPU_arithmetic_asian_option";
var GPU_basket_option_url = "/get_GPU_arithmetic_basket_option";

function get_GPU_arith_price(){
    $spot_price = $("#asian_spot_price").val();
    $volatility = $("#asian_volatility").val();
    $rate = $("#asian_rate").val();
    $maturity = $("#asian_maturity").val();
    $strike_price = $("#asian_strike_price").val();
    $observation_num = $("#asian_observation_num").val();
    $path_num = $("#asian_path_num").val();
    $control_variate = $("#asian_control_variate").val();
    $option_type = $("#asian_option_type").val();
    $quasi = $("#asian_quasi").prop('checked');
    url = GPU_arith_option_url;
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        observation_num:$observation_num,
        path_num:$path_num,
        control_variate:$control_variate,
        option_type:$option_type,
        quasi:$quasi};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            $output = dates.split(':');
            $("#asian_price").hide().fadeIn("slow");
            $("#asian_option_price").html($output[0]);
            $option_detail = "<br>Spot price : "+$spot_price+"; Volatility : "+$volatility
                    +";<br>"+"Risk-free interest Rate : "+$rate+"; Time to Maturity : "+$maturity
                    +";<br>"+"Strike price : "+$strike_price+"; Option Type : "+$option_type
                    +";<br>"+"Observation Number : "+$observation_num
                    +";<br>"+"Path Number : "+$path_num+"; Control Variate : "+$control_variate
                    +";<br> quasi : "+$quasi;

            //alert($option_detail);
            //$("#asian_option").html($option_detail);
            $("#asian_conf").hide().fadeIn("slow");
            $("#asian_option_conf").html($output[1]);
            //$("#get_asian_option").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
     });

}
function get_GPU_basket_arith_price(){
    $S1 = $("#basket_stock1_price").val();
    $S2 = $("#basket_stock2_price").val();
    $V1 = $("#basket_volatility1").val();
    $V2 = $("#basket_volatility2").val();
    $R = $("#basket_rate").val();
    $T = $("#basket_maturity").val();
    $K = $("#basket_strike_price").val();
    $rou = $("#basket_correlation").val();
    $path_num = $("#basket_path_num").val();
    $control_variate = $("#basket_control_variate").val();
    $option_type = $("#basket_option_type").val();
    $quasi = $("#basket_quasi").prop('checked');

    if (parseFloat($rou) > 1){
        alert("the abs of correlation must smaller than 1!");
        $("#basket_correlation").focus();
        return
    }
    url = GPU_basket_option_url;
    var data = {S1:$S1,
        S2:$S2,
        V1:$V1,
        V2:$V2,
        R:$R,
        T:$T,
        K:$K,
        rou:$rou,
        path_num:$path_num,
        control_variate:$control_variate,
        option_type:$option_type,
        quasi:$quasi};

    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $output = dates.split(':');
            $("#basket_price").hide().fadeIn("slow");
            $("#basket_option_price").html($output[0]);
            $option_detail = "Stock1 Price : "+$S1
            +";Stock2 Price : "+$S2+";<br>"+"Volatility1  : "+$V1
            +"; Volatility2 : "+$V2+";<br>"+"Risk-free interest Rate : "+$R
            +"; Time to Maturity : "+$T+";<br>"+"Strike price : "+$K
            +"; Correlation : "+$rou+";<br>"+"Path Number: "+$path_num
            +"; Control Variate: "+$control_variate+";<br>"
            +"Option Type : "+$option_type +"; quasi : "+$quasi;
            //$("#basket_option").html($option_detail);
            $("#basket_conf").hide().fadeIn("slow");
            $("#basket_option_conf").html($output[1]);
            //$("#get_basket_option").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
    });
}

function get_MC_euro_option() {
    $spot_price = $("#euro_spot_price").val();
    $volatility = $("#euro_volatility").val();
    $rate = $("#euro_rate").val();
    $maturity = $("#euro_maturity").val();
    $strike_price = $("#euro_strike_price").val();
    $option_type = $("#euro_option_type").val();
    $observation_num = $("#euro_observation_num").val();
    $path_num = $("#euro_path_num").val();
    url = MC_euro_option_url;
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        option_type:$option_type,
        observation_num:$observation_num,
        path_num:$path_num};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $output = dates.split(':');
            $("#euro_price").hide().fadeIn("slow");
            $("#euro_option_price").html($output[0]);
            $option_detail = "<br>Spot price : "+$spot_price+"; Volatility : "+$volatility
                    +";<br>"+"Risk-free interest Rate : "+$rate+"; Time to Maturity : "+$maturity
                    +";<br>"+"Strike price : "+$strike_price+"; Option Type : "+$option_type;
            //$("#euro_option").html($option_detail);
            //$("#bs_euro").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
    });
}

function get_GPU_euro_option() {
    $spot_price = $("#euro_spot_price").val();
    $volatility = $("#euro_volatility").val();
    $rate = $("#euro_rate").val();
    $maturity = $("#euro_maturity").val();
    $strike_price = $("#euro_strike_price").val();
    $option_type = $("#euro_option_type").val();
    $observation_num = $("#euro_observation_num").val();
    $path_num = $("#euro_path_num").val();
    $quasi = $("#euro_quasi").prop('checked');
    url = GPU_euro_option_url;
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        option_type:$option_type,
        observation_num:$observation_num,
        path_num:$path_num,
        quasi:$quasi};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $output = dates.split(':');
            $("#euro_price").hide().fadeIn("slow");
            $("#euro_option_price").html($output[0]);
            $option_detail = "<br>Spot price : "+$spot_price+"; Volatility : "+$volatility
                    +";<br>"+"Risk-free interest Rate : "+$rate+"; Time to Maturity : "+$maturity
                    +";<br>"+"Strike price : "+$strike_price+"; Option Type : "+$option_type;
            //$("#euro_option").html($option_detail);
            //$("#bs_euro").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
    });
}