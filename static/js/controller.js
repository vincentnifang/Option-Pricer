function get_geo_price(){
    $spot_price = $("#asian_spot_price").val();
    $volatility = $("#asian_volatility").val();
    $rate = $("#asian_rate").val();
    $maturity = $("#asian_maturity").val();
    $strike_price = $("#asian_strike_price").val();
    $observation_num = $("#asian_observation_num").val();
    $option_type = $("#asian_option_type").val();
    var url = "/get_geometric_asian_option";
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        observation_num:$observation_num,
        option_type:$option_type};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $("#asian_conf").hide();
            $("#asian_price").hide().fadeIn("slow");
            $("#asian_option_price").html(dates);
            $option_detail = "<br>Spot price : "+$spot_price+"; Volatility : "+$volatility
                    +";<br>"+"Risk-free interest Rate : "+$rate+"; Time to Maturity : "+$maturity
                    +";<br>"+"Strike price : "+$strike_price+"; Option Type : "+$option_type
                    +";<br>"+"Observation Number : "+$observation_num;
            //alert($option_detail)
            $("#asian_option").html($option_detail);
            //$("#get_asian_option").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
     });

}
function get_arith_price(){
    $spot_price = $("#asian_spot_price").val();
    $volatility = $("#asian_volatility").val();
    $rate = $("#asian_rate").val();
    $maturity = $("#asian_maturity").val();
    $strike_price = $("#asian_strike_price").val();
    $observation_num = $("#asian_observation_num").val();
    $path_num = $("#asian_path_num").val();
    $control_variate = $("#asian_control_variate").val();
    $option_type = $("#asian_option_type").val();
    var url = "/get_arithmetic_asian_option";
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        observation_num:$observation_num,
        path_num:$path_num,
        control_variate:$control_variate,
        option_type:$option_type};
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
                    +";<br>"+"Path Number : "+$path_num+"; Control Variate : "+$control_variate;

            //alert($option_detail);
            $("#asian_option").html($option_detail);
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

function get_basket_geo_price(){
    $S1 = $("#basket_stock1_price").val();
    $S2 = $("#basket_stock2_price").val();
    $V1 = $("#basket_volatility1").val();
    $V2 = $("#basket_volatility2").val();
    $R = $("#basket_rate").val();
    $T = $("#basket_maturity").val();
    $K = $("#basket_strike_price").val();
    $rou = $("#basket_correlation").val();
    $option_type = $("#basket_option_type").val();
    if (parseFloat($rou) > 1){
        alert("the abs of correlation must smaller than 1!");
        $("#basket_correlation").focus();
        return
    }
    var url = "/get_geometric_basket_option";
    var data = {S1:$S1,
        S2:$S2,
        V1:$V1,
        V2:$V2,
        R:$R,
        T:$T,
        K:$K,
        rou:$rou,
        option_type:$option_type};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $("#basket_conf").hide();
            $("#basket_price").hide().fadeIn("slow");
            $("#basket_option_price").html(dates);
            $option_detail = "Stock1 Price : "+$S1
            +";Stock2 Price : "+$S2+";<br>"+"Volatility1  : "+$V1
            +"; Volatility2 : "+$V2+";<br>"+"Risk-free interest Rate : "+$R
            +"; Time to Maturity : "+$T+";<br>"+"Strike price : "+$K
            +"; Correlation : "+$rou+";<br>"
            +"Option Type : "+$option_type;
            $("#basket_option").html($option_detail);
            //$("#get_basket_option").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//要刷新的div
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
    });
}

function get_basket_arith_price(){
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

    if (parseFloat($rou) > 1){
        alert("the abs of correlation must smaller than 1!");
        $("#basket_correlation").focus();
        return
    }
    var url = "/get_arithmetic_basket_option";
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
        option_type:$option_type};



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
            +"Option Type : "+$option_type;
            $("#basket_option").html($option_detail);
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
function get_bs_price(){
    $spot_price = $("#bs_spot_price").val();
    $volatility = $("#bs_volatility").val();
    $rate = $("#bs_rate").val();
    $maturity = $("#bs_maturity").val();
    $strike_price = $("#bs_strike_price").val();
    $option_type = $("#bs_option_type").val();
    var url = "/bs_euro";
    var data = {spot_price:$spot_price,
        volatility:$volatility,
        rate:$rate,
        maturity:$maturity,
        strike_price:$strike_price,
        option_type:$option_type};
    $.ajax({
        type : "get",
        //async : false,  //同步请求
        url : url,
        data : data,
        timeout:1000,
        success:function(dates){
            //alert(dates);
            $("#bs_price").hide().fadeIn("slow");
            $("#bs_option_price").html(dates);
            $option_detail = "<br>Spot price : "+$spot_price+"; Volatility : "+$volatility
                    +";<br>"+"Risk-free interest Rate : "+$rate+"; Time to Maturity : "+$maturity
                    +";<br>"+"Strike price : "+$strike_price+"; Option Type : "+$option_type;
            $("#bs_option").html($option_detail);
            //$("#bs_euro").closest('form').find("input[type=text], textarea").val("");
            //$("#out").html(dates);//
        },
        error: function() {
            alert("Pls filled correct numbers and try again!");
        }
    });
}