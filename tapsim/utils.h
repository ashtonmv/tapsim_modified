#ifndef TAPSIM_UTILS_H
#define TAPSIM_UTILS_H

/************************************************************************
*                                                                       *
* TAPSim - an atom probe data simulation program                        *
*                                                                       *
* Copyright (C) 2011 Christian Oberdorfer                               *
*                                                                       *
* This program is free software: you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation, either version 3 of the License, or any *
* any later version.                                                    *
*                                                                       *
* This program is distributed in the hope that it will be useful,       *
* but WITHOUT ANY WARRANTY; without even the implied warranty of        *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
* GNU General Public License for more details.                          *
*                                                                       *
* You should have received a copy of the GNU General Public License     *
* along with this program.  If not, see 'http://www.gnu.org/licenses'   *
*                                                                       *
************************************************************************/ 

#include <stdexcept>
#include <vector>
#include <cmath>

////////////////////////////////////////////////////////////////////////
//                                                                    //
// lu_decmp()                                                         //
//                                                                    //
// Source: Numerical Recipes in C, 2nd edition                        //
//                                                                    //
// Given a matrix lu[0..n-1][1..n-1], this routine replaces it by the //
// LU decomposition of a rowwise permutation of itself. lu and n are  //
// input. lu is output, arranged as with lower and upper tridiagonal  //
// elements joined in one matrix. By convention the digonale elements //
// of the  lower tridiagonal are set to  unity.  index[0..n-1] is an  //
// output  vector  that records the  row permutation effected  by the //
// partial pivoting;  d is  output as ± 1 depending  on whether  the  //
// number of row interchanges was even or odd, respectively.          //
// This routine  is used in  combination with  lu_solve() to  solve   //
// linear equations, invert a matrix or compute the determinant of a  //
// matrix.                                                            //
//                                                                    //
// *** Modified by Christian Oberdorfer                               //
// *** = returns false if the input matrix is singular                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

template<class REAL_TYPE, int n>
bool lu_decmp(REAL_TYPE *lu, int *index =0, REAL_TYPE *d =0)
{
	REAL_TYPE parity(1.0);
	REAL_TYPE ii[n];
	REAL_TYPE vv[n];

	for (int i = 0; i < n; i++)
	{
		REAL_TYPE big(0.0);
		for (int j = 0; j < n; j++)
		{
			const REAL_TYPE tmp = std::fabs(lu[i*n+j]);
			if (tmp > big) big = tmp;
		}
		
		if (big == 0.0) return false; // matrix is singular

		vv[i] = 1.0;
		vv[i] /= big;
	}

	int imax(0);

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < j; i++)
		{
			REAL_TYPE sum = lu[i*n+j];
			for (int k = 0; k < i; k++)
				sum -= lu[i*n+k] * lu[k*n+j];
			lu[i*n+j] = sum;
		}

		// ***

		REAL_TYPE big(0.0);
		for (int i = j; i < n; i++)
		{
			REAL_TYPE sum = lu[i*n+j];
			for (int k = 0; k < j; k++)
				sum -= lu[i*n+k] * lu[k*n+j];
			lu[i*n+j] = sum;

			// ***

			const REAL_TYPE dum = vv[i] * std::fabs(sum);

			if (dum >= big)
			{
				big = dum;
				imax = i;
			}
		}
		
		// ***

		if (j != imax)
		{
			for (int k = 0; k < n; k++)
				std::swap<REAL_TYPE>(lu[imax*n+k],lu[j*n+k]);

			parity *= -1.0;
			
			vv[imax] = vv[j];
		}

		ii[j] = imax;

		if (lu[j*n+j] == 0.0) return false; // matrix is singular

		if (j < n)
		{
			const REAL_TYPE dum = 1.0/lu[j*n+j];
			for (int i = j+1; i < n; i++)
				lu[i*n+j] *= dum;
		}
	}

	if (index != 0)
	{
		for (int i = 0; i < n; i++)
			index[i] = ii[i];
	}
	
	if (d != 0) *d = parity;
	
	return true;
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
// lu_solve()                                                         //
//                                                                    //
// Source: Numerical Recipes in C, 2nd edition                        //
//                                                                    //
// Solves the set of n linear equations A * X = B. Here LU is input,  //
// not as the matrix A but rather as its LU decomposition, determined //
// by the routine lu_decmp(). Index[0..n-1] is input as the permuta-  //
// tion vector returned by lu_decmp(). b[0..n-1] is input as the      //
// right-hand side vector B, and returns with the solution vector X.  //
// lu, n, and index are not modified by this routine and can be left  //
// in place for successive calls with different right-hand sides b.   //
// This routine takes into account the possibility that b will begin  //
// with many zero elements, so it is efficient for use in matrix      //
// inversion.                                                         //
//                                                                    //
// *** Modified by Christian Oberdorfer                               //
//                                                                    //
////////////////////////////////////////////////////////////////////////

template<class REAL_TYPE, int n>
void lu_solve(const REAL_TYPE* const lu, const int* const index, REAL_TYPE* b)
{
	int ii(0);

	for (int i = 0; i < n; i++)
	{
		const int ip = index[i];
		REAL_TYPE sum = b[ip];
		b[ip] = b[i];

		if (ii)
		{
			for (int j = ii-1; j < i; j++)
				sum -= lu[i*n+j] * b[j];
		}
		else
		{
			if (sum) ii = i + 1;
		}

		b[i] = sum;
	}
	
	for (int i = n-1; i >= 0; i--)
	{
		REAL_TYPE sum = b[i];
		for(int j = i+1; j < n; j++)
			sum -= lu[i*n+j] * b[j];

		b[i] = sum / lu[i*n+i];
	}
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
// Computation of the determinant from a square matrix                //
//                                                                    //
// Source: Numerical Recipes in C, 2nd edition                        //
//                                                                    //
// Parameters: lu,n and d is the output as obtained from lu_decmp()   //
//                                                                    //
// The determinant is the product of the diagonal elements in lu.     //
// The correct sign is given by the parity in d. - In order to avoid  //
// numeric underflow or overflow during the calculation, exponent and //
// mantisse are treated separately.                                   //
//                                                                    //
////////////////////////////////////////////////////////////////////////

template<class REAL_TYPE, int n>
REAL_TYPE lu_determinant(const REAL_TYPE* const lu, const REAL_TYPE d)
{
	REAL_TYPE x(d);
	int y(0);

	for (int i = 0; i < n; i++)
	{
		int tmp;
		x *= std::frexp(lu[i*n+i],&tmp);
		y += tmp;
	}

	return std::ldexp(x,y);
}

////////////////////////////////////////////////////////////////////////
//                                                                    //
// Reduction of a symmetric matrix to tridiagonal form by conducting  //
// the Householder reduction                                          //
//                                                                    //
// Source: Numerical Recipes in C, 2nd edition                        //
//                                                                    //
// Householder reduction of a real, symmetric matrix a[1...n][1...n]. //
// On output, a is replaced by the orthogonal matrix Q effecting the  //
// transformation. d[1..n] returns the diagonal elements of the tri-  //
// diagonal matrix, and e[1..n] the off-diagonal elements, with       //
// e[1]=0. Several statements, as noted in comments, can be ommitted  //
// if only eigenvalues are to be found, in which case a contains no   //
// useful information on output. Otherwise they are to be included.   //
//                                                                    //
////////////////////////////////////////////////////////////////////////

template<class REAL_TYPE, int n>
void householder(REAL_TYPE* a, REAL_TYPE* d, REAL_TYPE* e)
{
	int l,k,j,i;
	REAL_TYPE scale,hh,h,g,f;

	for (i=n-1;i>0;i--)
	{
		l=i-1;
		h=scale=0.0;
		if (l > 0)
		{
			for (k=0;k<l+1;k++)
				scale += std::fabs(a[i*n+k]);
			if (scale == 0.0)
				e[i]= a[i*n+l];
			else
			{
				for (k=0;k<l+1;k++)
				{
					a[i*n+k] /= scale;
					h += a[i*n+k]*a[i*n+k];
				}
				f=a[i*n+l];
				g=(f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i*n+l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++)
				{
					// Next statement can be omitted if eigenvectors not wanted
					a[j*n+i]=a[i*n+j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j*n+k]*a[i*n+k];
					for (k=j+1;k<l+1;k++)
						g += a[k*n+j]*a[i*n+k];
					e[j]=g/h;
					f += e[j]*a[i*n+j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++)
				{
					f=a[i*n+j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
				}
			}
		}
		else
			e[i]=a[i*n+l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	// wanted except for statement d[i]=a[i][i];

	for (i=0;i<n;i++)
	{
		l=i;

		if (d[i] != 0.0)
		{
			for (j=0;j<l;j++)
			{
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i*n+k]*a[k*n+j];
				for (k=0;k<l;k++)
					a[k*n+j] -= g*a[k*n+i];
			}
		}

		d[i]=a[i*n+i];
		a[i*n+i]=1.0;

		for (j=0;j<l;j++)
			a[j*n+i]=a[i*n+j]=0.0;
	}
}

/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Eigenvalues and eigenvectors of a tridiagonal matrix                //
//                                                                     //
// Source: Numerical Recipes in C, 2nd edition                         //
//                                                                     //
// QL algorithm with implicit shifts, to determine the eigenvalues     //
// and eigenvectors of a real, symmetric, tridiagonal matrix, or of a  //
// real, symmetric previously reduced by householder(). On input,      //
// d[0...n-1] contains the diagonal elements of the tridiagonal        //
// matrix. On output, it returns the eigenvalues. The vector e[0..n-1] //
// inputs the subdiagonal elements of the tridiagonal matrix, with     //
// e[0] arbitrary. On output e is destroyed. When finding only the     //
// eigenvalues, severel lines may be ommitted, as noted in the         //
// comments. If the eigenvectors of a tridiagonal matrix are desired,  //
// the matrix z[0...n-1][0...n-1] is input as the identity matrix. If  //
// the eigenvectors of a matrix that has been reduced by householder() //
// are required, then z is input as the matrix output by               //
// householder(). In either case, the kth column of z returns the      //
// normalized eigenvector corresponding to d[k].                       //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

template<class REAL_TYPE, int n>
void tri_eigen_qli(REAL_TYPE* d, REAL_TYPE* e, REAL_TYPE* z)
{
	int m,l,iter,i,k;
	REAL_TYPE s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++)
		e[i-1]=e[i];

	e[n-1]=0.0;

	for (l=0;l<n;l++)
	{
		iter=0;

		do
		{
			for (m=l;m<n-1;m++)
			{
				dd=std::fabs(d[m])+std::fabs(d[m+1]);
				if (std::fabs(e[m])+dd == dd) break;
			}

			if (m != l)
			{
				if (iter++ == 30) std::runtime_error("Too many iterations in tri_eigen_qli()");
     
				g=(d[l+1]-d[l])/(2.0*e[l]);

				{
					//r=pythag(g,1.0);

					const REAL_TYPE absa = std::fabs(g);
					const REAL_TYPE absb = std::fabs(1.0);

					REAL_TYPE pythag;

					if (absa > absb) 
						pythag =  absa*std::sqrt(1.0+(absb/absa)*(absb/absa));
					else
						pythag = (absb == 0.0 ? 0.0 : absb * std::sqrt(1.0+(absa/absb)*(absa/absb)));

					r = pythag;
				}

				//g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				// #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
				g = d[m]-d[l]+e[l]/(g+(g >= 0.0 ? std::fabs(r) : -std::fabs(r)));

				s = 1.0;
				c = 1.0;
				p = 0.0;

				for (i=m-1;i>=l;i--)
				{
					f=s*e[i];
					b=c*e[i];

					{
						// e[i+1] = (r=pythag(f,g));

						const REAL_TYPE absa = std::fabs(f);
						const REAL_TYPE absb = std::fabs(g);

						REAL_TYPE pythag;

						if (absa > absb) 
							pythag = absa*std::sqrt(1.0+(absb/absa)*(absb/absa));
						else
							pythag = ( absb == 0.0 ? 0.0 : absb * std::sqrt(1.0+(absa/absb)*(absa/absb)));

						e[i+1] = pythag;
						r = pythag;
					}

					if (r == 0.0)
					{
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i+1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					d[i+1] = g + (s * r);
					p = s * r;
					g = c * r - b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=0;k<n;k++)
					{
						f=z[k*n+i+1];
						z[k*n+i+1]=s*z[k*n+i]+c*f;
						z[k*n+i]=c*z[k*n+i]-s*f;
					}
				}

				if (r == 0.0 && i >= l) continue;
     
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}

#endif






