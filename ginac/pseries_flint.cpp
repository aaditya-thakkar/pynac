#include <flint/fmpq_polyxx.h>
#include <numeric>
#include <limits>
#include <stdexcept>
#include "ex.h"
#include "symbol.h"
#include "numeric.h"

namespace GiNaC {

using fp_t = flint::fmpq_polyxx;

/* Returns the sine series expansion of flint */
fp_t series_sine(const fp_t& s, const fp_t& var, unsigned int prec)
{
    return fp_t(sin_series(s,prec));                                      
}

using series_func_t = decltype(&series_sine);       /* type-defining series_func_t to the function type of series_sine() */

/* function map of the trig-functions */
static std::map<TypeID, series_func_t> function_map;
void build_function_map() {
    const symbol x = symbol("x");
    function_map[sin(x)->get_type_code()] = &series_sine;
} 

/* recursive function, which will convert the function argument in the series expansion and then pass it to the upper level function by recurrence */
fp_t _series(const ex &x, const symbol &var, unsigned int prec)
{
    if (is_exactly_a<numeric>(x)) {
        const numeric i1 = static_cast<const numeric &>(x);
        flint::fmpqxx i(i1.as_mpz().get_mpz_t())
        return fp_t(i);
    }

    if (is_exactly_a<symbol>(x)) {
        return pp_t(dynamic_cast<const symbol &>(x).get_name());
    }
    
    const fp_t var_p(var.get_name());
    if (is_exactly_a<function>(x)) {
        if (x.nops() > _ex1) {
            throw std::logic_error("Expansion of multivar functions not Implemented");
        }
        if (function_map.size() == 0) {
            build_function_map();
        } 
        series_func_t func = function_map[x.get_type_code()];
        if (func == nullptr) {
            throw std::logic_error(std::string("No expansion for this function: ")+ x.get_name());
        }
        auto series_inner = _series(x.op(0), var, prec);
        return func(std::move(series_inner), var_p, prec);
    }
}

}
