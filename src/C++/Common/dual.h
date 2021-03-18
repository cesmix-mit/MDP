// -*-c++-*-
//
// The MIT License (MIT)
// 
// Copyright (c) 2006 Jeffrey A. Fike
// Copyright (C) 2015 Michael Tesch tesch1 a gmail com
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//
// Version 0.0.2
//
// Thanks to Jeff Fike for publishing the dual original code upon
// which this file is based.
//
// The cxxdual package is available on github, bug reports and feature
// additions are welcome.
//
//  Docs:   http://tesch1.github.io/cxxduals
//  Source: https://github.com/tesch1/cxxduals
//
// #ifndef CXXDUALS_NO_COMPLEX
//
#ifndef LIB_CXXDUALS
#define LIB_CXXDUALS 1

#include <cmath>
#include <ctgmath>
#include <sstream>
#include <limits>
#include <type_traits>
#include <initializer_list>
#include <random>

#ifndef CXXDUALS_NO_COMPLEX
#include <complex>
#endif

// check if constexpr exists
#if __cplusplus >= 201103L || defined(_WIN32)
#define CXXDUALS_CONSTEXPR constexpr
#else
#define CXXDUALS_CONSTEXPR const
#endif

#ifdef __CUDACC__
// /////////////// Stuff missing in CUDA ///////////////////
#define DUAL_DEVICE_FUNC __host__ __device__
#define DUAL_STD_MATH(FUNC) using ::FUNC

#include <float.h>
#ifndef NPP_MIN_8U
#include <nppdefs.h>
#undef NV_NPPIDEFS_H
#endif
#include <math_constants.h>

template <typename T> struct numeric_limits;

template <> struct numeric_limits<float>
{
  __device__ __forceinline__ static float quiet_NaN() { return CUDART_NAN_F; };
  __device__ __forceinline__ static float infinity() { return CUDART_INF_F; };
  __device__ __forceinline__ static float epsilon() { return FLT_EPSILON; };
  __device__ __forceinline__ static float min() { return NPP_MINABS_32F; };
  __device__ __forceinline__ static float max() { return NPP_MAXABS_32F; };
};

template <> struct numeric_limits<double>
{
  __device__ __forceinline__ static double quiet_NaN() { return CUDART_NAN; };
  __device__ __forceinline__ static double infinity() { return CUDART_INF; };
  __device__ __forceinline__ static double epsilon() { return DBL_EPSILON; };
  __device__ __forceinline__ static double min() { return NPP_MINABS_64F; };
  __device__ __forceinline__ static double max() { return NPP_MAXABS_64F; };
};

template <> struct numeric_limits<short>
{
  __device__ __forceinline__ static short max() { return SHRT_MAX; };
};

template <typename _Tp>
bool isnormal(const _Tp & x)
{
  // maybe not exactly right?
  return !(x==0.0 || isnan(x) || isinf(x) || isfinite(x));
}

#else
// /////////////// Not CUDA ///////////////////

#define DUAL_DEVICE_FUNC
#define DUAL_STD_MATH(FUNC) using std::FUNC
#endif
// /////////////// End CUDA ///////////////////

/*!
 * dual-numbers implementation for calculation of differentials
 */
namespace cxxduals {

// forward declaration
template <typename _Tp> class dual;

// useful typedefs
#ifndef CXXDUALS_NO_TYPEDEFS
typedef dual<float> dualf;
typedef dual<double> duald;
typedef dual<long double> dualld;
#ifndef CXXDUALS_NO_COMPLEX
typedef dual<std::complex<float> > dualcf;
typedef dual<std::complex<double> > dualcd;
typedef dual<std::complex<long double> > dualcld;
#endif // CXXDUALS_NO_COMPLEX

#if __cplusplus >= 201103L || defined(_WIN32)
template <class _Tp> using hyperdual = dual<dual<_Tp> >;
typedef hyperdual<float> hyperdualf;
typedef hyperdual<double> hyperduald;
typedef hyperdual<long double> hyperdualld;
#ifndef CXXDUALS_NO_COMPLEX
typedef hyperdual<std::complex<float> > hyperdualcf;
typedef hyperdual<std::complex<double> > hyperdualcd;
typedef hyperdual<std::complex<long double> > hyperdualcld;
#endif // CXXDUALS_NO_COMPLEX
#endif // __cplusplus >= 201103L
#endif // CXXDUALS_NO_TYPEDEFS

} // cxxduals

//! muckin about in std:: to enable std::numeric_limits<> is allowed
//! by C++03 17.4.3.1/1, and C++11 18.3.2.3/1
#ifndef CXXDUALS_NO_LIMITS
namespace std {
#if 1
template<typename _Tp>
struct numeric_limits<cxxduals::dual<_Tp> > : public numeric_limits<_Tp> {
  typedef cxxduals::dual<_Tp> T;
  static CXXDUALS_CONSTEXPR bool is_specialized = true;
  static CXXDUALS_CONSTEXPR T min()           { return T(numeric_limits<_Tp>::min()); }
  static CXXDUALS_CONSTEXPR T max()           { return T(numeric_limits<_Tp>::max()); }
  static CXXDUALS_CONSTEXPR T epsilon()       { return T(numeric_limits<_Tp>::epsilon()); }
  static CXXDUALS_CONSTEXPR T round_error()   { return T(numeric_limits<_Tp>::round_error()); }
  static CXXDUALS_CONSTEXPR T infinity()      { return T(numeric_limits<_Tp>::infinity()); }
  static CXXDUALS_CONSTEXPR T quiet_NaN()     { return T(numeric_limits<_Tp>::quiet_NaN()); }
  static CXXDUALS_CONSTEXPR T signaling_NaN() { return T(numeric_limits<_Tp>::signaling_NaN()); }
  static CXXDUALS_CONSTEXPR T denorm_min()    { return T(numeric_limits<_Tp>::denorm_min()); }
};
#endif
// TODO: cv-specializations (C++11 18.3.2.3/2)

//! common_type specializations for dual
template <typename _Tp, typename _Up>
struct common_type<cxxduals::dual<_Tp>, cxxduals::dual<_Up> > {
  typedef cxxduals::dual<typename std::common_type<_Tp, _Up>::type> type;
};
} // std
#endif // CXXDUALS_NO_LIMITS

namespace cxxduals {

#ifndef CXXDUALS_NO_COMPLEX
/// check if something is complex
template <class T>
struct is_complex : std::false_type {};
template <class T>
struct is_complex<std::complex<T> > : public std::true_type {};
#endif

/// check if something is dual
template <class T>
struct is_dual : std::false_type {};
template <class T>
struct is_dual<dual<T> > : std::true_type {};

namespace internal {
/// Useful type trait extraction for templpated types, helper for _Tp
/// being a basic type (float, double, ...).  Used by dual_traits<>.
template <typename _Tp>
struct dual_traits_basic {
  typedef _Tp scalar_type;              ///< The intrinsic scalar type (double, float,...)
  typedef _Tp value_type;               ///< The type of part(), rpart() and epart()
  typedef _Tp basic_value_type;         ///< Used to differentiate complex/dual types
  static const int depth = 0;           ///< Depth of nested duals, a
                                        ///  dual<double> is 0, a
                                        ///  dual<dual<double>> is 1
  static const int num_elem = 1;        ///< How many scalar parts the dual has (2^depth)
  typedef std::false_type is_nested;    ///< Is this a nested dual?
  typedef std::false_type is_dual;      ///< Is _Tp even a dual?
};
}

/// Get trait information about a dual number (see
/// internal::dual_traits_basic for traits)
template <typename _Tp> struct dual_traits;

/// Specialization of dual_traits for float
template <> struct dual_traits<float> : public internal::dual_traits_basic<float> {};

/// Specialization of dual_traits for double
template <> struct dual_traits<double> : public internal::dual_traits_basic<double> {};

/// Specialization of dual_traits for long double
template <> struct dual_traits<long double> : public internal::dual_traits_basic<long double> {};

#ifndef CXXDUALS_NO_COMPLEX
/// Complex dual traits
template <typename _Tp>
struct dual_traits<std::complex<_Tp> > : public dual_traits<_Tp> {
  typedef typename dual_traits<_Tp>::scalar_type scalar_type;
  typedef _Tp value_type;
  typedef std::complex<_Tp> basic_value_type;
};
#endif // CXXDUALS_NO_COMPLEX

/// Get traits of duals
template <typename _Tp>
struct dual_traits<dual<_Tp> > {
  typedef typename dual_traits<_Tp>::scalar_type scalar_type;
  typedef _Tp value_type;
  typedef typename dual_traits<_Tp>::basic_value_type basic_value_type;
  static const int depth = dual_traits<_Tp>::depth + 1;
  static const int num_elem = dual_traits<_Tp>::num_elem * 2;
  typedef std::true_type is_dual;
  typedef typename dual_traits<_Tp>::is_dual is_nested;
};

namespace internal {
/// filter to exclude complex from templated arguments
template <typename _Tp> struct nocx_filter;
template <> struct nocx_filter<int> { typedef int test; };
template <> struct nocx_filter<float> { typedef int test; };
template <> struct nocx_filter<double> { typedef int test; };
template <> struct nocx_filter<long double> { typedef int test; };

/// Filter to enable tempated arg if a condition B is met
template< bool B, class T = void >
using enable_if_t = typename std::enable_if<B,T>::type;

/// Filter to only allow basic types in templated arg type
template <typename _Tp> struct arg_filter;
template <> struct arg_filter<int> { typedef int test; };
template <> struct arg_filter<float> { typedef int test; };
template <> struct arg_filter<double> { typedef int test; };
template <> struct arg_filter<long double> { typedef int test; };

template <typename _Tp> struct is_arithmetic : public std::is_arithmetic<_Tp> {
  typedef typename std::is_arithmetic<_Tp>::type type;
};

#ifndef CXXDUALS_NO_COMPLEX
#if 1
template <> struct arg_filter<std::complex<float> > { typedef int test; };
template <> struct arg_filter<std::complex<double> > { typedef int test; };
template <> struct arg_filter<std::complex<long double> > { typedef int test; };
#else
template <typename _Xp>
struct arg_filter<std::complex<_Xp> > : public arg_filter<_Xp> { typedef int test; };
#endif // 1
template <typename _Tp> struct is_arithmetic<std::complex<_Tp> > : public std::integral_constant<bool, true> {};
#endif // CXXDUALS_NO_COMPLEX

}

/// dual number class
template <typename _Tp>
class dual {

private:
  _Tp _f0, _f1;

public:

  /// Type of rpart() and epart()
  typedef _Tp value_type;

  // 
  typedef typename dual_traits<_Tp>::basic_value_type basic_value_type;
  typedef typename dual_traits<_Tp>::scalar_type scalar_type;
  static const int depth = dual_traits<dual<_Tp> >::depth;
  static const int num_elem = dual_traits<dual<_Tp> >::num_elem;

  /// Constructor - no initialization
  DUAL_DEVICE_FUNC
  CXXDUALS_CONSTEXPR dual()
    : _f0(), _f1() { }

  /// Constructor - specify real part, zero dual part
  DUAL_DEVICE_FUNC
  CXXDUALS_CONSTEXPR dual(const _Tp & f0)
    : _f0(f0), _f1(0) { }

  /// Constructor - specify real and dual parts
  DUAL_DEVICE_FUNC
  CXXDUALS_CONSTEXPR dual(const _Tp & f0, const _Tp & f1)
    : _f0(f0), _f1(f1) { }

  /// Construct from similar dual<> \todo restrict to depth<f0> == depth<_Tp>
  template <typename _Up,
            typename _Ux = typename internal::arg_filter<_Up>::test >
  DUAL_DEVICE_FUNC
  CXXDUALS_CONSTEXPR dual(const _Up & f0)
    : _f0(f0), _f1(0) { }

  /// Construct from similar dual<>
  /// \todo restrict to (depth<f0> == depth<_Tp> and depth<f1> == depth<_Tp>)
  template <typename _Up, typename _Vp, 
            typename _Ux = typename internal::arg_filter<_Up>::test,
            typename _Vx = typename internal::arg_filter<_Vp>::test >
  DUAL_DEVICE_FUNC
  CXXDUALS_CONSTEXPR dual(const _Up & f0, const _Vp & f1)
    : _f0(f0), _f1(f1) { }

#if __cplusplus >= 201103L || defined(_WIN32)
  /// List initializer, this is for allowing nested
  //  duals to have simpler initializers, ie for
  // dual<dual<double>> x{1,2,3,4}; // 1 + 2*e1 + 3*e2 + 4*e3
  DUAL_DEVICE_FUNC
  dual(std::initializer_list<basic_value_type> ll)
    : _f0(), _f1()
  {
#ifndef __CUDACC__
    if (ll.size() > num_elem)
      throw std::exception();
#endif
    int ii = 0;
    for (auto it = ll.begin(); it != ll.end(); it++, ii++)
      part(ii) = *it;
  }
#endif

  /// We're friends with all other duals (\todo ...of equal depth)
  template<typename _Up> friend class dual;

  /// Copy constructor from dual of diff type
  template<typename _Up,
           internal::enable_if_t<dual_traits<dual<_Up> >::depth <= depth>* = nullptr>
  DUAL_DEVICE_FUNC
  dual(const dual<_Up> & rhs)
    : _f0(rhs._f0), _f1(rhs._f1) { }

#if __cplusplus >= 201103L || defined(_WIN32)
  /// Real part
  DUAL_DEVICE_FUNC
  constexpr _Tp 
  rpart() const { return _f0; }

  /// Dual (epsilon) part
  DUAL_DEVICE_FUNC
  constexpr _Tp 
  epart() const { return _f1; }

private:
  /// Helper for `constexpr part(int) const`
  DUAL_DEVICE_FUNC
  constexpr basic_value_type
  private_part(int p, std::false_type) const {
    // p had better be either 0 or 1
    return p == 0 ? _f0 : _f1;
  }

  /// Helper for `constexpr part(int) const`
  DUAL_DEVICE_FUNC
  constexpr basic_value_type
  private_part(int p, std::true_type) const {
    return p < (num_elem / 2)
               ? _f0.part(p)
               : _f1.part(p - (dual_traits<dual<_Tp> >::num_elem / 2));
  }

public:
  /// Part-wise access, for nested duals
  DUAL_DEVICE_FUNC
  constexpr basic_value_type
  part(int p) const {
    return private_part(p, typename dual_traits<_Tp>::is_dual());
  }

#else

  /// Real part
  DUAL_DEVICE_FUNC
  inline const _Tp &
  rpart() const { return _f0; }

  /// Dual (epsilon) part
  DUAL_DEVICE_FUNC
  inline const _Tp &
  epart() const { return _f1; }

#endif

  /// Real part
  DUAL_DEVICE_FUNC
  inline _Tp &
  rpart() { return _f0; }

  /// Dual (epsilon) part
  DUAL_DEVICE_FUNC
  inline _Tp &
  epart() { return _f1; }

private:
  /// Helper for `inline part(int)`
  DUAL_DEVICE_FUNC
  inline basic_value_type &
  private_part(int p, std::false_type) {
    return p == 0 ? _f0 : _f1;
  }

  /// Helper for `inline part(int)`
  DUAL_DEVICE_FUNC
  inline basic_value_type &
  private_part(int p, std::true_type) {
    return p < (num_elem / 2)
               ? _f0.part(p)
               : _f1.part(p - (dual_traits<dual<_Tp> >::num_elem / 2));
  }

public:

  /// Part-wise access
  DUAL_DEVICE_FUNC
  inline basic_value_type &
  part(int p) {
    return private_part(p, typename dual_traits<_Tp>::is_dual());
  }

  /// Real part assignment
  DUAL_DEVICE_FUNC
  inline void
  rpart(_Tp f0) { _f0 = f0; }

  /// Dual (epsilon) part assignment
  DUAL_DEVICE_FUNC
  inline void
  epart(_Tp f1) { _f1 = f1; }

private:
  /// Helper for part-wise assignment
  DUAL_DEVICE_FUNC
  inline void
  private_part(int p, const basic_value_type & v, std::false_type) {
    if (p == 0)
      _f0 = v;
    else
      _f1 = v;
  }

  /// Helper for part-wise assignment
  DUAL_DEVICE_FUNC
  inline void
  private_part(int p, const basic_value_type & v, std::true_type) {
    if (p < (num_elem / 2))
      _f0.part(p, v);
    else
      _f1.part(p - (dual_traits<dual<_Tp> >::num_elem / 2), v);
  }

public:

  /// Set a single part value
  DUAL_DEVICE_FUNC
  inline void
  part(int p, const basic_value_type & v) {
    private_part(p, v, typename dual_traits<_Tp>::is_dual());
  }

  /// Assignment
  template <typename _Up,
            typename _Ux = typename internal::arg_filter<_Up>::test >
  DUAL_DEVICE_FUNC
  inline dual<_Tp> &
  operator=(const _Up & f0) { _f0 = f0; _f1 = _Tp(); return *this; }

  /// Assignment from another dual<> \todo restrict to depth<_Up> == depth<_Tp>
  template<typename _Up>
  DUAL_DEVICE_FUNC
  inline dual<_Tp> &
  operator=(const dual<_Up> & rhs) { _f0 = rhs._f0; _f1 = rhs._f1; return *this; }

  // Operations

  /// explicit cast to Scalar type -- same as rpart()
  DUAL_DEVICE_FUNC
  explicit operator _Tp () const { return rpart(); }

  /// unitary plus - passthrough
  DUAL_DEVICE_FUNC
  dual<_Tp>
  operator+() const { return *this; }

  /// unitary negation
  DUAL_DEVICE_FUNC
  dual<_Tp>
  operator-() const
  {
    return dual<_Tp>(-_f0, -_f1);
  }

  /// addition
  /**
   * \f[
   * (a + b \epsilon) + c
   * \f]
   */
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator+=(const _Tp & rhs)
  {
    _f0 += rhs;
    return *this;
  }

  /**
   * \f[
   * (a + b \epsilon) + (c + d \epsilon)
   * \f]
   */
  template<typename _Up>
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator+=(const dual<_Up> & rhs)
  {
    _f0 += rhs.rpart();
    _f1 += rhs.epart();
    return *this;
  }

  /// subtraction
  /**
   * \f[
   * (a + b \epsilon) - c
   * \f]
   */
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator-=(const _Tp & rhs)
  {
    _f0 -= rhs;
    return *this;
  }

  /**
   * \f[
   * (a + b \epsilon) - (c + d \epsilon)
   * \f]
   */
  template<typename _Up>
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator-=(const dual<_Up> & rhs)
  {
    _f0 -= rhs.rpart();
    _f1 -= rhs.epart();
    return *this;
  }

  /// multiplication
  /**
   * \f[
   * (a + b \epsilon) * c
   * \f]
   */
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator*=(const _Tp & rhs)
  {
    _f0 *= rhs;
    _f1 *= rhs;
    return *this;
  }

  /// multiplication
#ifndef CXXDUALS_NO_COMPLEX
  /**
   * Scalar multiplication when the value_type is complex
   * \f[
   * (a + (b+e*i) \epsilon) * c = (ac + (b+e*i)c \epsilon )
   * \f]
   */
  template<typename _Up>
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator*=(const typename std::complex<_Up> & rhs)
  {
    _f0 *= rhs;
    _f1 *= rhs;
    return *this;
  }
#endif // CXXDUALS_NO_COMPLEX

  /**
   * \f[
   * (a + b \epsilon) * (c + d \epsilon)
   * \f]
   *
   * \f[
   * 
   * \begin{pmatrix}
   * a & b \\
   * 0 & a \\
   * \end{pmatrix}
   *
   * \begin{pmatrix}
   * c & d \\
   * 0 & c \\
   * \end{pmatrix}
   * =
   * \begin{pmatrix}
   * ac & ad+bc \\
   * 0  & ac    \\
   * \end{pmatrix}
   * 
   * \f]
   */
  template<typename _Up>
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator*=(const dual<_Up> & rhs)
  {
    _Tp aa, bb;
    _Tp cc, dd;
    aa = _f0;
    bb = _f1;
    cc = rhs.rpart();
    dd = rhs.epart();
    _f0 = aa * cc;
    _f1 = aa * dd + bb * cc;
    return *this;
  }

  /// division
  /**
   * \f[
   * (a + b \epsilon) / c
   * \f]
   */
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator/=(const _Tp & rhs)
  {
    _f0 /= rhs;
    _f1 /= rhs;
    return *this;
  }

  /**
   * \f[
   * (a + b \epsilon) / (c + d \epsilon)
   * \f]
   *
   * \f[
   * 
   * \begin{pmatrix}
   * a & b \\
   * 0 & a \\
   * \end{pmatrix}
   *
   * \begin{pmatrix}
   * c & d \\
   * 0 & c \\
   * \end{pmatrix} ^ {-1}
   * =
   * \begin{pmatrix}
   * a & b \\
   * 0 & a \\
   * \end{pmatrix}
   * \frac{1}{c^2}
   * \begin{pmatrix}
   * c  & -d   \\
   * 0  &  c   \\
   * \end{pmatrix}
   * =
   * \frac{1}{c^2}
   * \begin{pmatrix}
   * ac & -ad + bc  \\
   * 0  & ac        \\
   * \end{pmatrix}
   * =
   * \begin{pmatrix}
   * a/c & \frac{bc - ad}{c^2}  \\
   * 0   & a/c                  \\
   * \end{pmatrix} \\
   *
   * = (a/c + \frac{bc - ad}{c^2} \epsilon)
   *
   * \f]

   */
  template<typename _Up>
  DUAL_DEVICE_FUNC
  dual<_Tp> &
  operator/=(const dual<_Up> & rhs)
  {
#if __cplusplus > 199711L || defined(_WIN32)
    typedef decltype(_Tp(1) * _Up(1)) higher_t;
#else
    typedef _Tp higher_t;
#endif
    //assert(rhs.rpart() != 0); // exclude non-field numbers
    dual<higher_t> tmp;
    tmp.rpart() = _f0 / rhs.rpart();
    tmp.epart() = (_f1 * rhs.rpart() - _f0 * rhs.epart()) / (rhs.rpart() * rhs.rpart());
    *this = tmp;
    return *this;
  }

};

template <typename _Up>
dual<_Up> copysign(dual<_Up> x, dual<_Up> y)
{
  x.rpart() = ::copysign(x.rpart(), y.rpart());
  return x;
}

/// Value extraction
//@{
/**
 * \f[
 * rpart(a + b \epsilon) = a
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline _Tp
rpart(const dual<_Tp> & d)
{
  return d.rpart();
}

/**
 * \f[
 * rpart(a) = a
 * \f]
 */
template <typename _Tp,
          typename _Tpf = typename internal::arg_filter<_Tp>::test>
DUAL_DEVICE_FUNC
inline _Tp
rpart(const _Tp & d)
{
  return d;
}

#ifndef CXXDUALS_NO_COMPLEX
/**
 * \f[
 * rpart((a + b \epsilon) + i*(c + d \epsilon)) = (a + c*i)
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline std::complex<_Tp>
rpart(const std::complex<dual<_Tp> > & d)
{
  return std::complex<_Tp>(d.real().rpart(), d.imag().rpart());
}
#endif

/**
 * \f[
 * epart(a + b \epsilon) = b
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline _Tp epart(const dual<_Tp> & d)
{
  return d.epart();
}

/**
 * \f[
 * epart(a) = 0
 * \f]
 */
template <typename _Tp,
          typename _Tpf = typename internal::arg_filter<_Tp>::test>
DUAL_DEVICE_FUNC
inline _Tp
epart(const _Tp & d)
{
  return _Tp(0);
}

//@}

/**
 * Random number in [a, b) - wrapper around std::default_random_engine
 * \f[
 * 
 * \f]
 */
template <typename _Tp,
          typename _Tpf = typename internal::nocx_filter<_Tp>::test>
DUAL_DEVICE_FUNC
void rand(_Tp & d, _Tp a = 0, _Tp b = 1)
{
  static std::default_random_engine generator;
  std::uniform_real_distribution<_Tp> distribution(a, b);
  d = distribution(generator);
}

#ifndef CXXDUALS_NO_COMPLEX
/**
 * Random complex number
 * \f[
 * 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
void rand(std::complex<_Tp> & d,
          typename dual_traits<std::complex<_Tp> >::scalar_type a = 0,
          typename dual_traits<std::complex<_Tp> >::scalar_type b = 1)
{
  _Tp x;
  rand(x, a, b);
  d.real(x);
  rand(x, a, b);
  d.imag(x);
}
#endif

/**
 * Random dual number with parts in [a, b)
 * \f[
 * 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
void rand(dual<_Tp> & d,
          typename dual_traits<dual<_Tp> >::scalar_type a = 0,
          typename dual_traits<dual<_Tp> >::scalar_type b = 1)
{
  typename dual_traits<dual<_Tp> >::scalar_type x;
  for (int ii = 0; ii < dual_traits<dual<_Tp> >::num_elem; ii++) {
    rand(x, a, b);
    d.part(ii, x);
  }
}

/**
 * Random dual number (used by Eigen::Random())
 * \f[
 * 
 * \f]
 */
template <typename _Tp,
          internal::enable_if_t<is_dual<_Tp>::value>* = nullptr>
DUAL_DEVICE_FUNC
_Tp
random(typename dual_traits<_Tp>::scalar_type a = 0,
       typename dual_traits<_Tp>::scalar_type b = 1)
{
  _Tp d;
  rand(d, a, b);
  return d;
}

// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };

// basic ops
#define DUALH_DEFINE_BASIC_OP_TEMPLATES(OP)                     \
  template <typename _Tp, typename _Up>                         \
  DUAL_DEVICE_FUNC                                              \
  inline typename std::common_type<dual<_Tp>, dual<_Up> >::type \
  operator OP (const dual<_Tp> & lhs, const dual<_Up> & rhs)    \
  {                                                             \
    typename std::common_type<dual<_Tp>, dual<_Up> >::type d = lhs;     \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline dual<_Tp>                                              \
  operator OP (const dual<_Tp> & lhs, const typename identity_t<_Tp>::type & rhs) \
  {                                                             \
    dual<_Tp> d = lhs;                                          \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline dual<_Tp>                                              \
  operator OP (const typename identity_t<_Tp>::type & lhs, const dual<_Tp> & rhs) \
  {                                                             \
    dual<_Tp> d = lhs;                                          \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline dual<_Tp>                                              \
  operator OP (const dual<_Tp> & lhs, const typename _Tp::value_type & rhs) \
  {                                                             \
    dual<_Tp> d = lhs;                                          \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline dual<_Tp>                                              \
  operator OP (const typename _Tp::value_type & lhs, const dual<_Tp> & rhs) \
  {                                                             \
    dual<_Tp> d = lhs;                                          \
    d OP##= rhs;                                                \
    return d;                                                   \
  }

#ifndef CXXDUALS_NO_COMPLEX
#define DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(OP)                  \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline std::complex<dual<_Tp>>                                     \
  operator OP (const std::complex<dual<_Tp>> & lhs, const std::complex<_Tp> & rhs) \
  {                                                             \
    std::complex<dual<_Tp>> d = lhs;                            \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline std::complex<dual<_Tp>>                                        \
  operator OP (const std::complex<_Tp> & lhs, const std::complex<dual<_Tp> > & rhs) \
  {                                                             \
    std::complex<dual<_Tp>> d = lhs;                            \
    d OP##= rhs;                                                \
    return d;                                                   \
  }\
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline std::complex<dual<_Tp>>                                     \
  operator OP (const std::complex<dual<_Tp>> & lhs, const _Tp & rhs) \
  {                                                             \
    std::complex<dual<_Tp>> d = lhs;                            \
    d OP##= rhs;                                                \
    return d;                                                   \
  }                                                             \
  template <typename _Tp>                                       \
  DUAL_DEVICE_FUNC                                              \
  inline std::complex<dual<_Tp>>                                      \
  operator OP (const _Tp & lhs, const std::complex<dual<_Tp> > & rhs) \
  {                                                             \
    std::complex<dual<_Tp>> d = lhs;                            \
    d OP##= rhs;                                                \
    return d;                                                   \
  }
#else
#define DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(OP)
#endif

//@{
///  Addition
DUALH_DEFINE_BASIC_OP_TEMPLATES(+)
DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(+)
//@}

//@{
///  Subtraction
DUALH_DEFINE_BASIC_OP_TEMPLATES(-)
DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(-)
//@}

//@{
///  Multiplication
DUALH_DEFINE_BASIC_OP_TEMPLATES(*)
DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(*)
//@}

//@{
///  Division
DUALH_DEFINE_BASIC_OP_TEMPLATES(/)
DUALH_DEFINE_BASIC_OP_TEMPLATES_CX(/)
//@}

#undef DUALH_DEFINE_BASIC_OP_TEMPLATES

//
#if 0
// fix emacs auto formatting below this point for me
  ;
#endif

/**
 * \f[
 * (a + b \epsilon) ^ c = a^c + b c a^{c-1} \epsilon
 * \f]
 *
 * by the function derivative argument:
 * 
 * \f[
 * f(x) = x ^ c \\
 * f'(x) = c x ^ {c-1} \\
 * x = a + b \epsilon \\
 * f(a + b \epsilon) = f(a) + b f'(a) \epsilon \\
 *                   = a^c + b c a ^ {c-1} \epsilon \\
 * \f]
 */
template <typename _Tp, typename _Up,
          internal::enable_if_t<internal::is_arithmetic<_Up>::value>* = nullptr>
DUAL_DEVICE_FUNC
dual<_Tp>
pow(const dual<_Tp> & xx, const _Up & cc)
{
  DUAL_STD_MATH(pow);
  DUAL_STD_MATH(abs);
  DUAL_STD_MATH(exp);
#if 0
  dual<_Tp> temp;
  _Tp deriv, xval, tol;
  xval = xx.rpart();
  // TODO- should use numeric traits of _Tp instead of 1e-15
  tol = _Tp(1e-15);
  if (abs(xval) > 0 && abs(xval) < abs(tol)) {
    xval = xx.rpart() / (abs(xx.rpart()) / tol);
    //if (xval >= 0)
    //  xval = tol;
    //if (xval < 0)
    //  xval = -tol;
  }
  deriv = cc * pow(xval, (cc - _Tp(1.0)));
  temp.rpart() = pow(xx.rpart(), cc);  //Use actual x value, only use tol for derivs
  temp.epart() = xx.epart() * deriv;
#endif
  dual<_Tp> temp(pow(xx.rpart(), cc),
                 xx.epart() * _Tp(cc) * _Tp(pow(xx.rpart(), cc - _Up(1))));
  //temp.rpart() = pow(xx.rpart(), cc);
  //temp.epart() = xx.epart() * _Tp(cc) * _Tp(pow(xx.rpart(), cc - _Up(1)));
  return temp;
}

#if 1
/**
 * \f[
 * a ^ {c + d \epsilon} = a^c + d * a ^ c * \log(a) \epsilon
 * \f]
 *
 * by the function derivative argument:
 * 
 * \f[
 * f(y) = a ^ y \\
 * f'(y) = a ^ y \log (a) \\
 * y = c + d \epsilon \\
 * f(c + d \epsilon) = f(c) + d f'(c) \epsilon \\
 *                   = a^c + d a^c \log(a) \epsilon \\
 *
 * \f]
 */
template <typename _Tp, typename _Up,
          internal::enable_if_t<internal::is_arithmetic<_Tp>::value>* = nullptr>
DUAL_DEVICE_FUNC
dual<_Up>
pow(const _Tp & aa, const dual<_Up> & yy)
{
  DUAL_STD_MATH(pow);
  DUAL_STD_MATH(log);
  dual<_Up> temp;
  temp.rpart() = pow(aa, yy.rpart());
  temp.epart() = yy.epart() * log(aa) * pow(aa, yy.rpart());
  return temp;
}

/**
 * \f[
 * (a + b \epsilon)^{(c + d \epsilon)} = a^c + (b c a^{c-1} + d \log(a) a^c)\epsilon
 * \f]
 *
 * by the function derivative argument:
 * 
 * \f[
 * f(x,y) = x ^ y \\
 * f'(x,y) = 
 * x = a + b \epsilon \\
 * y = c + d \epsilon \\
 * f(x,y) = \\
 *
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
pow(const dual<_Tp> & xx, const dual<_Tp> & yy)
{
  DUAL_STD_MATH(pow);
  DUAL_STD_MATH(log);
  dual<_Tp> temp;
  temp.rpart() = pow(xx.rpart(), yy.rpart());
  temp.epart() =
    xx.epart() * yy.rpart() * pow(xx.rpart(), yy.rpart() - _Tp(1)) +
    yy.epart() * log(xx.rpart()) * pow(xx.rpart(), yy.rpart());
  return temp;
}
#endif

/**
 * \f[
 * \exp(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
exp(const dual<_Tp> & x)
{
  DUAL_STD_MATH(exp);
  return dual<_Tp>(exp(x.rpart()),
                   exp(x.rpart()) * x.epart());
}

/**
 * \f[
 * \log(a + b \epsilon) =  
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
log(dual<_Tp> x)
{
  DUAL_STD_MATH(log);
  _Tp deriv1;
  deriv1 = x.epart() / x.rpart();
  return dual<_Tp>(log(x.rpart()), deriv1);
}

/**
 * \f[
 * \sin(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline dual<_Tp>
sin(dual<_Tp> x)
{
  DUAL_STD_MATH(sin);
  DUAL_STD_MATH(cos);
  dual<_Tp> temp;
  _Tp funval, deriv;
  funval = sin(x.rpart());
  deriv = cos(x.rpart());
  temp.rpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

/**
 * \f[
 * \cos(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline dual<_Tp>
cos(dual<_Tp> x)
{
  DUAL_STD_MATH(sin);
  DUAL_STD_MATH(cos);
  dual<_Tp> temp;
  _Tp funval, deriv;
  funval = cos(x.rpart());
  deriv = -sin(x.rpart());
  temp.rpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

/**
 * \f[
 * \tan(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
inline dual<_Tp>
tan(dual<_Tp> x)
{
  DUAL_STD_MATH(tan);
  dual<_Tp> temp;
  _Tp funval, deriv;
  funval = tan(x.rpart());
  deriv  = funval*funval + 1.0;
  temp.rpart() = funval;
  temp.epart() = deriv*x.epart();
  return temp;
}

/**
 * \f[
 * \mathrm{asin}(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
asin(dual<_Tp> x)
{
  DUAL_STD_MATH(asin);
  DUAL_STD_MATH(sqrt);
  dual<_Tp> temp;
  _Tp funval, deriv1, deriv;
  funval = asin(x.rpart());
  deriv1 = 1.0 - x.rpart()*x.rpart();
  deriv = 1.0 / sqrt(deriv1);
  temp.rpart() = funval;
  temp.epart() = deriv*x.epart();
  return temp;
}

/**
 * \f[
 * \mathrm{acos}(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
acos(dual<_Tp> x)
{
  DUAL_STD_MATH(acos);
  DUAL_STD_MATH(sqrt);
  dual<_Tp> temp;
  _Tp funval, deriv1, deriv;
  funval = acos(x.rpart());
  deriv1 = 1.0 - x.rpart() * x.rpart();
  deriv = 1.0 / sqrt(deriv1);
  temp.rpart() = funval;
  temp.epart() = -deriv*x.epart();
  return temp;
}

/**
 * \f[
 * \mathrm{atan}(a + b \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
atan(dual<_Tp> x)
{
  DUAL_STD_MATH(atan);
  dual<_Tp> temp;
  _Tp funval, deriv1, deriv;
  funval = atan(x.rpart());
  deriv1 = 1.0 + x.rpart() * x.rpart();
  deriv = 1.0 / deriv1;
  temp.rpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

/**
 * \f[
 * \mathrm{atan2}(a + b \epsilon, c + d \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
atan2(dual<_Tp> y, dual<_Tp> x)
{
  DUAL_STD_MATH(atan2);
  dual<_Tp> temp;
  _Tp funval, deriv1, deriv;
  funval = atan2(y.rpart(), x.rpart());
  // unsure from here on...
  deriv1 = 1.0 + x.rpart() * x.rpart();
  deriv = 1.0 / deriv1;
  temp.rpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

/**
 * \f[
 * \sqrt{a + b \epsilon} =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
sqrt(dual<_Tp> x)
{
  DUAL_STD_MATH(pow);
  return pow(x, 0.5);
  //return pow(x, (typename dual<_Tp>::scalar_type)0.5);
}

/**
 * \f[
 * \max(a + b \epsilon, c + d \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
max(dual<_Tp> x1, dual<_Tp> x2)
{
  return x1.rpart() >= x2.rpart() ? x1 : x2;
}

/**
 * \f[
 * \max(a + b \epsilon, c) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
max(dual<_Tp> x1, _Tp x2)
{
  return x1.rpart() >= x2 ? x1 : dual<_Tp>(x2);
}

/**
 * \f[
 * \max(a, c + d \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
max(_Tp x1, dual<_Tp> x2)
{
  return x1 >= x2.rpart() ? dual<_Tp>(x1) : x2;
}

/**
 * \f[
 * \min(a + b \epsilon, c + d \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
min(dual<_Tp> x1, dual<_Tp> x2)
{
  return x1.rpart() <= x2.rpart() ? x1 : x2;
}

/**
 * \f[
 * \min(a + b \epsilon, c) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
min(dual<_Tp> x1, _Tp x2)
{
  return x1.rpart() <= x2 ? x1 : dual<_Tp>(x2);
}

/**
 * \f[
 * \min(a, c + d \epsilon) =
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
min(_Tp x1, dual<_Tp> x2)
{
  return x1 <= x2.rpart() ? dual<_Tp>(x1) : x2;
}

#ifndef CXXDUALS_NO_COMPLEX
/**
 * Complex Conjugation
 * \f[
 * 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<std::complex<_Tp> >
conj(const dual<std::complex<_Tp> > & x)
{
  // there is actually no derivative of conjugate()
  DUAL_STD_MATH(conj);
  return dual<std::complex<_Tp>>(conj(x.rpart()), conj(x.epart()));
}
#endif

/**
 * Dual Conjugation
 * \f[
 * i = \sqrt{-1} \\
 * \mathrm{real}( (e+f*i) + (g+h*i) \epsilon ) = e + g \epsilon
 * \f]
 *
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
dconj(const dual<_Tp> & x)
{
  return dual<_Tp>(x.rpart(), -x.epart());
}

#ifndef CXXDUALS_NO_COMPLEX
/**
 * Coupled Conjugation
 * \f[
 * 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<std::complex<_Tp> >
cconj(const dual<std::complex<_Tp> > & x)
{
  DUAL_STD_MATH(conj);
  return dual<std::complex<_Tp>>(conj(x.rpart()), -conj(x.epart()));
}

/**
 * Dual-Complex Conjugation
 * \f[
 * 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<std::complex<_Tp> >
dcconj(const dual<std::complex<_Tp> > & x)
{
  DUAL_STD_MATH(conj);
  return conj(x.rpart()) * (1 - x.rpart() / x.epart());
}

/**
 * Anti-Dual Conjugation
 * \f[
 * adconj(
 * \f]
 *
 * Messelmi, Farid. “DUAL-COMPLEX NUMBERS AND THEIR HOLOMORPHIC FUNCTIONS,” n.d., 12.
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<std::complex<_Tp> >
adconj(const dual<std::complex<_Tp> > & x)
{
  return dual<std::complex<_Tp> >(x.epart(), -x.rpart());
}
#endif

/**
 * \f[
 * i = \sqrt{-1} \\
 * \mathrm{real}( (e+f*i) + (g+h*i) \epsilon ) = e + g \epsilon
 * \f]
 *
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<typename _Tp::value_type>
real(const dual<_Tp> & x)
{
  // todo - dont just make things up
  return dual<typename _Tp::value_type>(real(x.rpart()),
                                        real(x.epart()));
}

/**
 * \f[
 * i = \sqrt{-1} \\
 * \mathrm{imag}( (e+f*i) + (g+h*i) \epsilon ) = f + h \epsilon
 * \f]
 *
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<typename _Tp::value_type>
imag(const dual<_Tp> & x)
{
  // todo - dont just make things up
  return dual<typename _Tp::value_type>(imag(x.rpart()),
                                        imag(x.epart()));
}

/**
 * \f[
 * |a + b \epsilon| = 
 * \f]
 *
 * by the xfunction derivative argument: (for real-valued functions)
 * 
 * \f[
 * f(x) = \sqrt{(\mathrm{real}(x)^2 + \mathrm{imag}(x)^2)} \\
 * f'(x) = \frac{x}{|x|} f'(x) \\
 * x = a + b \epsilon \\
 * f(a + b \epsilon) = f(a) + b f'(a) \epsilon \\
 *
 * \f]
 *
 * and for complex-valued functions:
 * \f[
 * ?
 * \f]
 *
 */
template <typename _Tp,
          typename _Tpf = typename internal::nocx_filter<_Tp>::test>
DUAL_DEVICE_FUNC
inline
dual<_Tp>
abs(const dual<_Tp> & x)
{
  DUAL_STD_MATH(abs);
  // according to ...
  return sqrt(dconj(x) * x);
  
  // how ceres does it
  // return abs(x.rpart()) == x.rpart() ? x : -x;

  // by way of arguing f' and
  // https://math.stackexchange.com/questions/1235365/find-the-derivative-of-absolute-value-using-the-chain-rule
  //return dual<_Tp>(abs(x.rpart()),
  //                 (x.rpart() / abs(x.rpart())) * x.epart() * x.epart());
}

#ifndef CXXDUALS_NO_COMPLEX
template <typename _Tp>
DUAL_DEVICE_FUNC
inline
std::complex<_Tp>
abs(const dual<std::complex<_Tp> > & x)
{
  DUAL_STD_MATH(conj);
  // return abs(x.rpart()) == x.rpart() ? x : -x;
  return rpart(conj(x) * x);
  //return x.rpart();
  //return dual<typename dual<_Tp>::scalar_type>(abs(x.rpart()),
  //                                             (x.rpart() / abs(x.rpart())) * x.epart() * x.epart());
}
#endif

/**
 * \f[
 * |a + b \epsilon|^2 = x
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
abs2(const dual<_Tp> & x)
{
  return x * x;
}

/**
 * \f[
 * \mathrm{arg}(a + b \epsilon) = 
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
_Tp
arg(const dual<_Tp> & x)
{
  // x.rpart() != 0
  return x.epart() / x.rpart();
}

/**
 * \f[
 * \mathrm{norm}(a + b \epsilon) = | a + b \epsilon | = a
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
_Tp
norm(const dual<_Tp> & x)
{
  return x.rpart();
}


/**
 * \f[
 *      
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
ceil(const dual<_Tp> & x)
{
  DUAL_STD_MATH(ceil);
  DUAL_STD_MATH(numeric_limits);
  _Tp c = ceil(x.rpart());
  return dual<_Tp>(c,
                   c == x.rpart() ?
                   numeric_limits<_Tp>::infinity()
                   : _Tp(0));
}

/**
 * \f[
 *      
 * \f]
 */
template <typename _Tp>
DUAL_DEVICE_FUNC
dual<_Tp>
floor(const dual<_Tp> & x)
{
  DUAL_STD_MATH(floor);
  DUAL_STD_MATH(numeric_limits);
  _Tp f = floor(x.rpart());
  return dual<_Tp>(f,
                   f == x.rpart()
                   ? numeric_limits<_Tp>::infinity()
                   : _Tp(0));
}

// fp classification
/// Are any of the parts NaN?
template <typename _Tp>
DUAL_DEVICE_FUNC
bool isnan(const dual<_Tp> & x)
{
  DUAL_STD_MATH(isnan);
  bool isit = true;
  for (int ii = 0; isit && ii < dual<_Tp>::num_elem; ii++)
    isit = isit && isnan(x.part(ii));
  return isit;
}

/// Are all of the parts not Inf?
template <typename _Tp>
DUAL_DEVICE_FUNC
bool isfinite(const dual<_Tp> & x)
{
  DUAL_STD_MATH(isfinite);
  bool isit = true;
  for (int ii = 0; isit && ii < dual<_Tp>::num_elem; ii++)
    isit = isit && isfinite(x.part(ii));
  return isit;
}

/// Are any of the parts Inf?
template <typename _Tp>
DUAL_DEVICE_FUNC
bool isinf(const dual<_Tp> & x)
{
  DUAL_STD_MATH(isinf);
  bool isit = true;
  for (int ii = 0; isit && ii < dual<_Tp>::num_elem; ii++)
    isit = isit && isinf(x.part(ii));
  return isit;
}

/// Are all of the parts isnormal()?
template <typename _Tp>
DUAL_DEVICE_FUNC
bool isnormal(const dual<_Tp> & x)
{
  DUAL_STD_MATH(isnormal);
  bool isit = true;
  for (int ii = 0; isit && ii < dual<_Tp>::num_elem; ii++)
    isit = isit && isnormal(x.part(ii));
  return isit;
}

/// comparison
#define DUALH_COMPARISON_OP(OP)                                 \
  template <typename _Tp, typename _Up>                         \
  DUAL_DEVICE_FUNC                                              \
  inline bool                                                   \
  operator OP (const dual<_Tp> & lhs, const dual<_Up> & rhs)    \
  {                                                             \
    return lhs.rpart() OP rhs.rpart();                                  \
  }                                                                     \
  template <typename _Tp, typename _Up,                                 \
            typename std::enable_if<internal::is_arithmetic<_Up>{},     \
                                    int>::type = 0>                     \
  DUAL_DEVICE_FUNC                                                  \
  inline bool                                                       \
  operator OP (const dual<_Tp> & lhs, const _Up & rhs)              \
  {                                                                 \
    return lhs.rpart() OP rhs;                                          \
  }                                                                     \
  template <typename _Tp, typename _Up,                                 \
            typename std::enable_if<internal::is_arithmetic<_Up>{},     \
                                    int>::type = 0>                     \
  DUAL_DEVICE_FUNC                                              \
  inline bool                                                   \
  operator OP (const _Up & lhs, const dual<_Tp> & rhs)          \
  {                                                             \
    return lhs OP rhs.rpart();                                  \
  }

DUALH_COMPARISON_OP(>)
DUALH_COMPARISON_OP(<)
DUALH_COMPARISON_OP(==)
DUALH_COMPARISON_OP(>=)
DUALH_COMPARISON_OP(<=)
DUALH_COMPARISON_OP(!=)

#undef DUALH_COMPARISON_OP
#if 0
// fix emacs auto formatting below this point for me
  ;
#endif
#undef DUAL_DEVICE_FUNC

//#include <iomanip>
template<typename _Tp, typename _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits> &
operator<<(std::basic_ostream<_CharT, _Traits> & os, const cxxduals::dual<_Tp> & rhs)
{
  using namespace cxxduals;
#if 0
  // print nested duals
  std::basic_ostringstream<_CharT, _Traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << "(" << rhs.rpart()
    << " + e" << dual_traits<dual<_Tp> >::depth << "*" << rhs.epart()
    << ")";
  return os << s.str();
#elif 1
  // print dual parts by epsilon-grade
  std::basic_ostringstream<_CharT, _Traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << "(";
  for (int p = 0; p < rhs.num_elem; p++) {
    if (p)
      s << " + e" << p << "*";
    s << rhs.part(p);
  }
  s << ")";
  return os << s.str();
#else
  // print as matrix (wrong)
  std::basic_ostringstream<_CharT, _Traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << "(";
  int side = rhs.num_elem;
  for (int i = 0; i < side; i++) {
    for (int j = 0; j < i; j++)
      //s << "         " << 0 << " ";
      s << std::setw(10) << 0 << " ";
    for (int j = i; j < side; j++) {
      s << std::setw(10) << rhs.part(j) << " ";
    }
    if (i < side-1)
      s << "\n ";
  }
  s << ")";
  return os << s.str();
#endif
}

} // namespace cxxduals

//////////////////////////////////////////////////////////////////////
// Eigen Support
//
#ifndef EIGEN_VERSION_AT_LEAST
#define EIGEN_VERSION_AT_LEAST(...) 0
#endif
#if defined(CXXDUALS_EIGEN) || EIGEN_VERSION_AT_LEAST(3, 3, 0)

namespace Eigen {

// This allows using cxxduals::dual in Eigen
template<typename _Scalar>
struct NumTraits<cxxduals::dual<_Scalar> > : GenericNumTraits<_Scalar>
{
  typedef cxxduals::dual<typename NumTraits<_Scalar>::Real> Real;
  //typedef cxxduals::dual<_Scalar> Real;
  //typedef _Scalar Real;
  typedef cxxduals::dual<typename NumTraits<_Scalar>::NonInteger> NonInteger;
  typedef cxxduals::dual<_Scalar> Nested;

  enum {
    IsInteger           =   NumTraits<_Scalar>::IsInteger,
    IsSigned            =   NumTraits<_Scalar>::IsSigned,
    IsComplex           =   0,
    //IsComplex           =   NumTraits<_Scalar>::IsComplex,
    RequireInitialization = NumTraits<_Scalar>::RequireInitialization,
    ReadCost            = 2 * NumTraits<_Scalar>::ReadCost,
    AddCost             = 2 * NumTraits<_Scalar>::AddCost,
    MulCost             = 4 * NumTraits<_Scalar>::MulCost + 2 * NumTraits<_Scalar>::AddCost
  };

  EIGEN_DEVICE_FUNC
  static inline Real epsilon()          { return Real(NumTraits<_Scalar>::epsilon()); }
  EIGEN_DEVICE_FUNC
  static inline Real dummy_precision()  { return Real(NumTraits<_Scalar>::dummy_precision()); }
};

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T, cxxduals::dual<T>, BinaryOp > {
  typedef cxxduals::dual<T> ReturnType;
};

#ifndef CXXDUALS_NO_COMPLEX
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T, cxxduals::dual<std::complex<T> >, BinaryOp > {
  typedef cxxduals::dual<std::complex<T> > ReturnType;
};
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T, std::complex<cxxduals::dual<T> >, BinaryOp> {
  typedef std::complex<cxxduals::dual<T> > ReturnType;
};
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<T>, std::complex<cxxduals::dual<T> >, BinaryOp > {
  typedef std::complex<cxxduals::dual<T> > ReturnType;
};
#endif

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<cxxduals::dual<T>, T, BinaryOp> {
  typedef cxxduals::dual<T> ReturnType;
};

#ifndef CXXDUALS_NO_COMPLEX
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<cxxduals::dual<std::complex<T> >, T, BinaryOp> {
  typedef cxxduals::dual<std::complex<T> > ReturnType;
};
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<cxxduals::dual<T> >, T, BinaryOp> {
  typedef std::complex<cxxduals::dual<T> > ReturnType;
};
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<cxxduals::dual<T> >, std::complex<T>, BinaryOp> {
  typedef std::complex<cxxduals::dual<T> > ReturnType;
};
#endif

} // namespace Eigen

namespace cxxduals {

/** define a custom template unary functor
 * use it like this: m2 = m1.unaryExpr(CwiseRpartOp<double>());
 */
template<typename Scalar>
struct CwiseRpartOp {
  const Scalar operator()(const Scalar & x) const { return x; }
  const Scalar operator()(const cxxduals::dual<Scalar> & x) const { return rpart(x); }
#ifndef CXXDUALS_NO_COMPLEX
  const std::complex<Scalar> operator()(const std::complex<Scalar> & x) const { return x; }
  const std::complex<Scalar> operator()(const std::complex<cxxduals::dual<Scalar> > & x) const {
    return std::complex<Scalar>(rpart(x.real()), rpart(x.imag()));
  }
#endif
};

/** define a custom template unary functor
 * use it like this: m2 = m1.unaryExpr(CwiseEpartOp<double>());
 */
template<typename Scalar>
struct CwiseEpartOp {
  const Scalar operator()(const cxxduals::dual<Scalar> & x) const { return epart(x); }
#ifndef CXXDUALS_NO_COMPLEX
  const std::complex<Scalar> operator()(const std::complex<Scalar> & x) const { return x; }
  const std::complex<Scalar> operator()(const std::complex<cxxduals::dual<Scalar> > & x) const {
    return std::complex<Scalar>(epart(x.real()), epart(x.imag()));
  }
#endif
};

}

#endif // CXXDUALS_EIGEN

// End Eigen Support
//////////////////////////////////////////////////////////////////////

#ifndef CXXDUALS_NO_COMPLEX

namespace cxxduals {

/// Make working with std::complex<> nubmers suck less... allow promotion.
#define COMPLEX_OPS(OP)                                                 \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }                                                                     \
                                                                        \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }

COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#if 0
  ;
#endif

}

#endif // CXXDUALS_NO_COMPLEX

#endif // LIB_CXXDUALS