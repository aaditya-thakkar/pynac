#include "ex.h"
#include "symbol.h"
#include "numeric.h"
#include "function.h"
#include "power.h"
#include <string>
#include <flint/fmpq_polyxx.h>

namespace GiNaC {

using fp_t = flint::fmpq_polyxx;

/* Returns the sine series expansion of flint */
fp_t series_sin(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(sin_series(s,prec));                                      
}

fp_t series_cos(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(cos_series(s,prec));                                      
}

fp_t series_tan(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(tan_series(s,prec));                                      
}

fp_t series_log(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(log_series(s,prec));                                      
}

fp_t series_exp(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(exp_series(s,prec));                                      
}

fp_t series_asin(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(asin_series(s,prec));                                      
}

fp_t series_atan(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(atan_series(s,prec));                                      
}

fp_t series_sinh(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(sinh_series(s,prec));                                      
}

fp_t series_cosh(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(cosh_series(s,prec));                                      
}

fp_t series_tanh(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(tanh_series(s,prec));                                      
}

fp_t series_asinh(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(asinh_series(s,prec));                                      
}

fp_t series_atanh(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(atanh_series(s,prec));                                      
}

using series_func_t = decltype(&series_sin);       /* type-defining series_func_t to the function type of series_sine() */

/* function map of the trig-functions */
static std::map<std::string, series_func_t> function_map;
void build_function_map() {
    function_map["sin"] = &series_sin;
    function_map["cos"] = &series_cos;
    function_map["tan"] = &series_tan;
    function_map["log"] = &series_log;
    function_map["exp"] = &series_exp;
    function_map["asin"] = &series_asin;
    function_map["atan"] = &series_atan;
    function_map["sinh"] = &series_sinh;
    function_map["cosh"] = &series_cosh;
    function_map["tanh"] = &series_tanh;
    function_map["asinh"] = &series_asinh;
    function_map["atanh"] = &series_atanh;
} 

static int _helper(const ex &x, const symbol &var, unsigned int prec)
{
    if (function_map.size() == 0) {
        build_function_map();
    }
    /* In case of numeric or symbolic argument, returns 1 */
    if (is_exactly_a<numeric>(x) || is_exactly_a<symbol>(x)) {
        return 1;
    }
    /* in case of function arguments if multivar case or upmapped function is encountered, 0 and -1 are returned */
    if (is_exactly_a<function>(x)) {
        if(x.nops() > 1) {
            return 0;
        }
        function f = ex_to<function>(x);
        series_func_t func = function_map[f.get_name()];
        if (func == nullptr) {
            return -1;
        }
    }
    /* calling _helper for the next internal level */
    const int rec = _helper(x.op(0), var, prec);
    /* returns 0 and -1 if anytime in the recurrence 0 or -1 have been returned, and 1 is returned if anytime in the recursion 1 has been returned */
    if (rec == -1) {
        return -1;
    }
    else if (rec == 0) {
        return 0;
    }
    else {
        return 1;
    }
} 

/* recursive function, which will convert the function argument in the series expansion and then pass it to the upper level function by recurrence */
fp_t _series(const ex &x, const symbol &var, unsigned int prec)
{
    const int h = _helper(x, var, prec);

    if (h == 0) {
        throw std::runtime_error("Expansion of multivar functions not Implemented");
    }
    if (h == -1) {
        throw std::runtime_error(std::string("No expansion for this function: "));
    }
    if (is_exactly_a<numeric>(x)) {
        const numeric i1 = ex_to<numeric>(x);
        std::string ss("1  ");
        std::string ss1("1  ");
        if (i1.is_rational()) {
           const numeric n = i1.numer();
           const numeric d = i1.denom();
           int num = n.to_int();
           int den = d.to_int();
           std::string strn = std::to_string(num);
           std::string strd = std::to_string(den);
           ss.append(strn);
           ss1.append(strd);
           const char *c = ss.c_str();
           const char *c1 = ss1.c_str();
           fp_t c2(fp_t(c)/fp_t(c1)); 
           return fp_t(c2);   
        }
        
        int a = i1.to_int();
        std::string str = std::to_string(a);
        ss.append(str);
        const char *c = ss.c_str();
        return fp_t(c);
    }

    if (is_exactly_a<symbol>(x)) {
        return fp_t("2  0 1");
    }
    
    const fp_t var_p("2  0 1");
    if (is_exactly_a<function>(x)) {
        series_func_t func = function_map[ex_to<function>(x).get_name()];
        auto series_inner = _series(x.op(0), var, prec);
        return func(std::move(series_inner), var_p, prec);
    }

    if (is_exactly_a<power>(x)) {
        const power p = ex_to<power>(x);
        const ex bs = p.get_basis();
        const ex expn = p.get_exponent();
        if (is_exactly_a<numeric>(expn)) {
            const numeric e = ex_to<numeric>(expn);
            if (e.is_rational()) {
                const numeric n = e.numer();
                const numeric d = e.denom();
                int num = n.to_int();
                int den = d.to_int();
                const fp_t ps_root(series_nthroot(_series(bs, var, prec), den, var_p, prec));   //series_nthroot() --> coming soon! :P
                if (num > 0) {
                    return fp_t(ps_root.pow(unsigned(num)));
                }
                else if (num < 0) {
                    return fp_t(ps_root.inv_series(prec).pow(unsigned(-num)));
                }
                return fp_t("1  0");
            } 
            int expon = e.to_int();
            const fp_t pbasis(_series(bs, var, prec));
            if (expon > 0) {
                return fp_t(pbasis.pow(unsigned(expon)));
            }
            else if (expon < 0) {
                return fp_t(pbasis.inv_series(prec).pow(unsigned(-expon)));
            }
            return fp_t("1  0");  
        }
        
    }
return fp_t("1  0");   
}

}
