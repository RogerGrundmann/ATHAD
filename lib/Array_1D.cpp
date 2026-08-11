/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 1D arrays
*/

#include <cassert>
#include <cmath>
#include <iostream>

#include "Array_1D.h"

using namespace std;

#define MAXSIZE 361

// create an Array_1D with specified size and initial value
Array_1D::Array_1D(int n, double val) : z(NULL){
    initArray_1D(n, val);
}
/*
*
*/
Array_1D::~Array_1D(){
    delete []z;
}
/*
*
*/
void Array_1D::initArray_1D(int mm, double cc){
    if(!z){//when z is null
        assert(mm <= MAXSIZE);
        this->mm = mm;
        z = new double[mm];
    }else{
        assert(mm == this->mm);
    }
    for(int i = 0; i < mm; i++){
        z[i] = cc;
    }
}
/*
*
*/
void Array_1D::Coordinates(int mm, double z0, double dz){
    assert(mm == this->mm); // FIXME: just until we remove mm throughout
    z[0] = z0;
    for(int l = 1; l < mm; l++){
        z[l] = z[l - 1] + dz;
    }
}
/*
*
*/
// Stretched coordinates — largest step at i=0, smallest at i=mm-1.
// Formula:  z[i] = z0 + L · (1 − sinh(β·(1 − i/(mm−1))) / sinh(β))
// where L = (mm−1)·dz.  At i=0: z=z0;  at i=mm-1: z=z0+L.
// Larger β increases the stretching ratio (surface step / bottom step = cosh(β)).
void Array_1D::StretchedCoordinates(int mm, double z0, double dz, double beta){
    assert(mm == this->mm);
    const double L        = (mm - 1) * dz;
    const double sinh_b   = sinh(beta);
    for (int i = 0; i < mm; i++) {
        const double xi = (double)i / (double)(mm - 1);
        z[i] = z0 + L * (1.0 - sinh(beta * (1.0 - xi)) / sinh_b);
    }
}
/*
*
*/
void Array_1D::printArray_1D(int mm){
    assert(mm == this->mm); // FIXME: just until we remove mm throughout
    cout.precision(6);
    cout.setf(ios::fixed);
    cout << endl;
    if(mm == 41) cout << " rad = i-direction in " << mm-1 << "  radial steps "  << endl;
    if(mm == 181) cout << " theta = j-direction in " << mm-1 << "  latitudinal steps "  << endl;
    if(mm == 361) cout << " phi = k-direction in " << mm-1 << "  longitudinal steps "  << endl;
    cout << endl;
    for(int i = 0; i < mm; i++){
        cout.width(6);
        cout.fill(' ');
        cout << z[i] << " ";
    }
    cout << endl;
}
