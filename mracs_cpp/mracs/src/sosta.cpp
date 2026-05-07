// Second Order STAtistics (SOSTA)of MRACS code, which include:
// (1) two-point correlation (with or without volume-average) function
// (2) density variance with top-hat filter
#include"mracs.h"

#define R0 0.5           // Mpc/h
#define R1 50.           // Mpc/h
#define NUMTEST 10

int main()
{
    read_parameter();
    std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    std::vector<Particle> p, p0; for(Galaxy i : g) p.push_back({i.x, i.y, i.z, 1.});
    std::vector<double> r_log, xi_r, var_r;
    
    for(int i = 0; i < NUMTEST; ++i) {
        r_log.push_back(R0 * pow((R1/R0), static_cast<double>(i)/NUMTEST));
    }
    auto s   = sfc(p);
    auto sc  = sfc_r2c(s);
    
    force_kernel_type(0);
    for(int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        xi_r.push_back(inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1);
        delete[] c;
    }
    force_kernel_type(1);   
    for( int i = 0; i < NUMTEST; ++i){
        auto w = wfc(r_log[i], 0);
        auto c = convol_c2r(sc, w);
        var_r.push_back(inner_product(c,c,GridNum)/pow(p.size()*4./3*M_PI*pow(r_log[i]/SimBoxL,3),2)/GridNum-1);
        delete[] c;
    }
    for(auto x : r_log) std::cout << x << ", "; std::cout << std::endl;
    for(auto x : xi_r)  std::cout << x << ", "; std::cout << std::endl;
    for(auto x : var_r) std::cout << sqrt(x) << ", "; std::cout << std::endl;
}
