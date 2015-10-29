#pragma once
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <eigen3/Eigen/Dense>

#ifdef DEBUG_BUILD
#   define DEBUG_CERR(x) do { std::cerr << x; } while (0)
#else
#   define DEBUG_CERR(x)
#endif

#ifdef PRECISION_LONG_DOUBLE
    typedef long double real;
#elif PRECISION_FLOAT
    typedef float real;
#else
    typedef double real;
#endif

typedef unsigned uint;
typedef std::complex<real> complex;
typedef Eigen::Matrix<real,-1,-1>    MatrixXr;
typedef Eigen::Matrix<complex,-1,-1> MatrixXc;
typedef Eigen::Matrix<real,-1,1>     VectorXr;
typedef Eigen::Matrix<complex,-1,1>  VectorXc;

typedef struct {
    MatrixXr A;
    MatrixXr widths;
    real norm;
} CGaussian;

typedef std::vector<CGaussian> Basis;

constexpr real pi    = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
constexpr real twopi = 2*pi;
constexpr real hbar  = 1;

//Wrapper to make a function return to a pointer supplied as an argument
template<typename F, typename T, typename... P>
void void_wrapper(F *func, T* out, P... params)
{
    *out = (*func)(params...);
}

//Same thing for member functions
template<class C, typename F, typename T, typename... P>
void void_member_wrapper(F C::* func, C* c, T* out, P... params)
{
    *out = (c->*func)(params...);
}

//Take a function returning non-void and run it as a thread, with
//output directed through a pointer
template<typename F, typename T, typename... P>
std::thread threadify(F *func, T* out, P... params)
{
    return std::thread(void_wrapper<F,T,P...>,func,out,params...);

}

template<class C, typename F, typename T, typename... P>
std::thread threadify_member(F C::* func, C* c, T* out, P... params)
{
    return std::thread(void_member_wrapper<C,F,T,P...>,func,c,out,params...);

}

namespace std {

    template<>
    struct hash<pair<string,string> > {
        size_t operator()(pair<string,string> const& ps) const {
            const size_t a ( hash<string>()(ps.first)  );
            const size_t b ( hash<string>()(ps.second) );
            return a^b;
        }
    };

    template<>
    struct equal_to<pair<string,string> > {
        bool operator()(pair<string,string> const& p1, pair<string,string> const& p2) const {
            if ( (p1.first == p2.first && p1.second == p2.second) ||
                 (p1.first == p2.second && p1.second == p2.first) ) {
                return true;
            } else {
                return false;
            }
        }
    };

}

#endif
