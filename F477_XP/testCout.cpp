

#include <iostream>

extern "C"
{

    void __mgeoidun_MOD_f477( double *tElapsed, double VecFlat[], double VecFlon[], double VecU[], int *num, int *nOrder );
    void __mgeoidun_MOD_testf477( double *tElapsed, double VecFlat[], double VecFlon[], double VecU[], int *num, int *nOrder, double resVec360[], double resVecIntpt[] );
    void __mgeoidun_MOD_init_const1( int *nOrder );
    void __intptfunc_MOD_intpt( double VecFlat[], double VecFlon[], double VecU[], int *num );

}

using namespace std;

int testIntpt()
{
    double flat[3], flon[3], vecu[3] = {0.0};

    //     data VecFlat /0.0, 90.0, -90.0/
    //   data VecFlon / 0.0, 180.0, -180.0/
    flat[0] = 0.0;
    flat[1] = 89.0;
    flat[2] = -89.0;

    flon[0] = 0.0;
    flon[1] = -179.0;
    flon[2] = -359.0;
    int num = 3;

    __intptfunc_MOD_intpt(flat, flon, vecu, &num);

    std::cout << "vecu : "<< vecu << std::endl;

    return 0;
}


int testF477(){

    double flat[3], flon[3], vecu[3] = {0.0},
     dElapsedT, res360[3]= {0.0}, resIntpt[3]= {0.0};

    //     data VecFlat /0.0, 90.0, -90.0/
    //   data VecFlon / 0.0, 180.0, -180.0/
    flat[0] = 0.0;
    flat[1] = 89.0;
    flat[2] = -89.0;

    flon[0] = 0.0;
    flon[1] = 179.0;
    flon[2] = 359.0;
    int num = 3;
    int nOrder = 30;


    __mgeoidun_MOD_init_const1( &nOrder );
    __mgeoidun_MOD_f477(&dElapsedT, flat, flon, vecu, &num, &nOrder);

    std::cout << "vecu : " << vecu[1] <<
     " dElapseT =  " << dElapsedT<< std::endl;


    __mgeoidun_MOD_testf477(&dElapsedT, flat, flon, vecu, &num, &nOrder, res360, resIntpt);


    std::cout << "vecu : " << vecu[1] <<
     " dElapseT =  " << dElapsedT<< std::endl;

    return 0;
}

int main(){

    testF477();
    testIntpt();

    return 0;
}



