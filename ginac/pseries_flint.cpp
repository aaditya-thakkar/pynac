#include <flint/fmpq_polyxx.h>

namespace GiNaC {

using fp_t = flint::fmpq_polyxx;

fp_t series_sine(const fp_t& s, const fp_t& var, unsigned int prec)
{
    fp_t p1;
    p1.set_zero();

    if(s == var) {
        // case: sin(x)
        flint::fmpqxx prod(0);
        for (unsigned int i=0; i<prec/2; i++) {
            const short j = 2*i + 1;
            if (i != 0) {
                prod *= 1-j;    
            }
            prod *= j;
            p1 += (1/prod) * var.pow(j);           //pow()=>will be implemented
        }
        return p1;
    }
    
    return fp_t(sin_series(s,prec));
    

}
