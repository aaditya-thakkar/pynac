#include <flint/fmpq_polyxx.h>
#include <numeric>
#include <limits>
#include <stdexcept>
#include "ex.h"
#include "symbol.h"
#include "numeric.h"
#include "function.h"

namespace GiNaC {

using fp_t = flint::fmpq_polyxx;

/* Returns the sine series expansion of flint */
fp_t series_sine(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(sin_series(s,prec));                                      
}

using series_func_t = decltype(&series_sine);       /* type-defining series_func_t to the function type of series_sine() */

/* function map of the trig-functions */
static std::map<function, series_func_t> function_map;
void build_function_map() {
    function_map[sin] = &series_sin;
   /* function_map[cos] = &series_cos;
    function_map[tan] = &series_tan;
    function_map[cot] = &series_cot;
    function_map[sec] = &series_sec;
    function_map[csc] = &series_csc;
    function_map[log] = &series_log;
    function_map[exp] = &series_exp;
    function_map[asin] = &series_asin;
    function_map[acos] = &series_acos;
    function_map[atan] = &series_atan;
    function_map[acot] = &series_acot;
    function_map[asec] = &series_asec;
    function_map[acsc] = &series_acsc;
    function_map[sinh] = &series_sinh;
    function_map[cosh] = &series_cosh;
    function_map[tanh] = &series_tanh;
    function_map[coth] = &series_coth;
    function_map[sech] = &series_sech;
    function_map[csch] = &series_csch;
    function_map[asinh] = &series_asinh;
    function_map[acosh] = &series_acosh;
    function_map[atanh] = &series_atanh;
    function_map[acoth] = &series_acoth;
    function_map[asech] = &series_asech;
    function_map[acsch] = &series_acsch;*/
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
        series_func_t func = function_map[f];
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
    if (h == 1) {
        if (is_exactly_a<numeric>(x)) {
            const numeric i1 = ex_to<numeric>(x);
            flint::fmpqxx i(i1);
            return fp_t(i);
        }

        if (is_exactly_a<symbol>(x)) {
            return fp_t(ex_to<symbol>(x).get_name());
        }
    
        const fp_t var_p(var.get_name());
        if (is_exactly_a<function>(x)) {
            series_func_t func = function_map[ex_to<function>(x)];
            auto series_inner = _series(x.op(0), var, prec);
            return func(std::move(series_inner), var_p, prec);
        }
    }
    else if (h == 0) {
        throw std::runtime_error("Expansion of multivar functions not Implemented");
    }
    else if (h == -1) {
        throw std::runtime_error(std::string("No expansion for this function: "));
    }
}


}
