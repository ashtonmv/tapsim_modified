#ifndef PREDICATES_H
#define PREDICATES_H

#ifdef __cplusplus

#include <stdexcept>
#include <string>

namespace Predicates
{
	extern "C"
	{
		void exactinit();

		  float orient2d_single(const float*, const float*, const float*);

		  // sign is according to the left-hand-rule (>0 if points are in clock-wise order)
		  float orient3d_single(const float*, const float*, const float*, const float*);

		  float incircle_single(const float*, const float*, const float*, const float*);
		  float insphere_single(const float*, const float*, const float*, const float*, const float*);

		  double orient2d_double(const double*, const double*, const double*);

		  // sign is according to the left-hand-rule (>0 if points are in clock-wise order)
		  double orient3d_double(const double*, const double*, const double*, const double*);

		  double incircle_double(const double*, const double*, const double*, const double*);
		  double insphere_double(const double*, const double*, const double*, const double*, const double*);
	}

	// ***

	template<class REAL_TYPE>
	inline REAL_TYPE orient2d(const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*)
	{
		throw std::runtime_error("orient2d(): template parameters must be 'float' or 'double'!");
	}

	template<class REAL_TYPE>
	inline REAL_TYPE orient3d(const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*)
	{
		throw std::runtime_error("orient3d(): template parameters must be 'float' or 'double'!");
	}

	template<class REAL_TYPE>
	inline REAL_TYPE incircle(const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*)
	{
		throw std::runtime_error("incircle(): template parameters must be 'float' or 'double'!");
	}

	template<class REAL_TYPE>
	inline REAL_TYPE insphere(const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*, const REAL_TYPE*)
	{
		throw std::runtime_error("insphere(): template parameters must be 'float' or 'double'!");
	}

	// ***

	template<>
	inline float orient2d<float>(const float* a, const float* b, const float* c)
	{
		return orient2d_single(a,b,c);
	}

	template <>
	inline float orient3d<float>(const float* a, const float* b, const float* c, const float* d)
	{
		return orient3d_single(a,b,c,d);
	}

	template<>
	inline float incircle<float>(const float* a, const float* b, const float* c, const float* d)
	{
		return incircle_single(a,b,c,d);
	}

	template<>
	inline float insphere<float>(const float* a, const float* b, const float* c, const float* d, const float* e)
	{
		return insphere_single(a,b,c,d,e);
	}

	// ***

	template<>
	inline double orient2d<double>(const double* a, const double* b, const double* c)
	{
		return orient2d_double(a,b,c);
	}

	template <>
	inline double orient3d<double>(const double* a, const double* b, const double* c, const double* d)
	{
		return orient3d_double(a,b,c,d);
	}

	template<>
	inline double incircle<double>(const double* a, const double* b, const double* c, const double* d)
	{
		return incircle_double(a,b,c,d);
	}

	template<>
	inline double insphere<double>(const double* a, const double* b, const double* c, const double* d, const double* e)
	{
		return insphere_double(a,b,c,d,e);
	}
}
#endif

#ifndef __cplusplus
	void exactinit();

	float orient2d_single(float*, float*, float*);

	// sign is according to the left-hand-rule (>0 if points are in clock-wise order)
	float orient3d_single(float*, float*, float*, float*);

	float incircle_single(float*, float*, float*, float*);
	float insphere_single(float*, float*, float*, float*, float*);

	double orient2d_double(double*, double*, double*);

	// sign is according to the left-hand-rule (>0 if points are in clock-wise order)
	double orient3d_double(double*, double*, double*, double*);

	double incircle_double(double*, double*, double*, double*);
	double insphere_double(double*, double*, double*, double*, double*);
#endif

#endif