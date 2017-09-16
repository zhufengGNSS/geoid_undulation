
#include <iostream>

//extern "C"{
//    void intptFunc_intpt( double *VecFlat, double *VecFlon, double *VecU );
//}

extern "C" void __intptfunc_MOD_intpt( double VecFlat[], double VecFlon[], double VecU[], int *num );

using namespace std;
int testIntpt2()
{
    double flat[3], flon[3], vecu[3] = {0.0};

    //     data VecFlat /0.0, 90.0, -90.0/
    //   data VecFlon / 0.0, 180.0, -180.0/
    flat[0] = 0.0;
    flat[1] = 89.0;
    flat[2] = -89.0;

    flon[0] = 0.0;
    flon[1] = 179.0;
    flon[2] = 359.0;
    int num = 3;

    __intptfunc_MOD_intpt(flat, flon, vecu, &num);

    std::cout << "vecu : "<< vecu << std::endl;

    return 0;
}
