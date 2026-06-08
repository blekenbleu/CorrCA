/* Matrix and vector classes.
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>
    (C) 2014 Victoria Rudakova <vicrucann@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef ARRAY_H // Do nothing if not included from array.h

namespace libNumerics {

/// Constructor
template <typename T>
vector<T>::vector(int m)
: matrix<T>(m, 1)
{}

template <typename T>
vector<T>::vector(void)
: matrix<T>()
{}

/// 1-vector constructor.
//template <typename T>
//vector<T>::vector(T x)
//: matrix<T>(1,1)
//{
//    this->p[0] = x;
//}

/// 2-vector constructor.
template <typename T>
vector<T>::vector(T x, T y)
: matrix<T>(2,1)
{
    this->p[0] = x;
    this->p[1] = y;
}

/// 3-vector constructor.
template <typename T>
vector<T>::vector(T x, T y, T z)
: matrix<T>(3,1)
{
    this->p[0] = x;
    this->p[1] = y;
    this->p[2] = z;
}

/// Copy constructor
template <typename T>
vector<T>::vector(const vector<T>& v)
: matrix<T>(v)
{}

/// initialize existing Vector
template <typename T>
void vector<T>::init(int m)
{
    alloc(m_rows = m, m_cols = 1);
}

template <typename T>
void vector<T>::init(T value, int m)
{
    alloc(m_rows = m, m_cols = 1);
	T *x = p + m;
    for(T *i = p; i < x; i++)
        *i = value;
}

/* Vector of indices
vector<int> vector<int>::index(int m)
{
    vector<int> V(m);
    for(int i = V.size()-1; i >= 0; i--)
        V.p[i] = i;
    return V;
} */

/// Assignment operator
template <typename T>
vector<T> &vector<T>::operator=(const vector<T> &v)
{
    matrix<T>::operator=(v);
    return *this;
}

/// Access the \a i-th coefficient.
template <typename T>
inline T vector<T>::operator[] (int i) const
{
    assert(i >= 0 && i < this->nElements());
    return this->p[i];
}

/// Access the \a i-th coefficient.
template <typename T>
inline T& vector<T>::operator[] (int i)
{
    assert(i >= 0 && i < this->nElements());
    return this->p[i];
}

/// Multiply a vector by scalar.
/// \param a a scalar.
template <typename T>
vector<T> vector<T>::operator*(T a) const
{
    vector<T> v(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        v.p[i] = a*this->p[i];
    return v;
}

/// Divide a vector by scalar.
/// \param a a scalar.
template <typename T>
inline vector<T> vector<T>::operator/(T a) const
{
    return operator*( (T)1/a );
}

/// Addition of vectors.
template <typename T>
vector<T> vector<T>::operator+(const vector<T>& v) const
{
    assert(this->m_rows == v.m_rows);
    vector<T> sum(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        sum.p[i] = this->p[i] + v.p[i];
    return sum;
}

/// Addition of constant value to each element of the vector.
template <typename T>
vector<T> vector<T>::operator+(T a) const
{
    vector<T> sum(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        sum.p[i] = this->p[i] + a;
    return sum;
}

/// Subtraction of vectors.
template <typename T>
vector<T> vector<T>::operator-(const vector<T> &v) const
{
    assert(this->m_rows == v.m_rows);
    vector<T> dif(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        dif.p[i] = this->p[i] - v.p[i];
    return dif;
}

/// Subtraction a constant value from each element of the vector.
template <typename T>
vector<T> vector<T>::operator-(T a) const
{
    vector<T> dif(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        dif.p[i] = this->p[i] - a;
    return dif;
}

/// Opposite of vector.
template <typename T>
vector<T> vector<T>::operator-() const
{
    vector<T> v(this->m_rows);
    for(int i = this->m_rows-1; i >= 0; i--)
        v.p[i] = -this->p[i];
    return v;
}

/// Vector times matrix.
template <typename T>
matrix<T> vector<T>::operator*(const matrix<T>& m) const
{
    return matrix<T>::operator*(m);
}

/// Dot product operator.
template <typename T>
T vector<T>::operator*(const vector<T>& v) const
{
    assert(v.size() == this->size());
#if (_WIN32 && (!_WIN64))
	return ddot_asm(this->p, v.p, this->size()); 
#else
	T res = (T)0;
    for(int i = v.nrow()-1; i >= 0; i--)
      res += this->p[i] * v[i];
    return res;
#endif
}

/// Diagonal matrix defined by its diagonal vector.
template <typename T>
matrix<T> vector<T>::diag() const
{
    matrix<T> d(this->m_rows, this->m_rows);
    d = (T)0;
    for(int i = this->m_rows-1; i >= 0; i--)
        d(i,i) = this->p[i];
    return d;
}

/// Square L^2 norm of vector.
template <typename T>
T vector<T>::qnorm() const
{
    T q = (T)0;
    for(int i = this->m_rows-1; i >= 0; i--)
        q += this->p[i]*this->p[i];
	return q;
}

/// SubvectorRef from \a i0 to \a i1.
template <typename T>
vectorRef<T> vector<T>::copyRef(int i0, int i1) const
{
	assert(0 <= i0 && i0 <= i1 && i1 <= this->m_rows);
	return vectorRef<T>(i1-i0+1, this->p+i0);
}

/// Subvector without row \a i0.
template<typename T>
void vector<T>::without(int i0, const vector<T> &v)
{
	int rows = v.m_rows;
	if (i0 > v.m_rows)
        i0 = v.m_rows;
    else {
		if (0 > i0)
        	i0 = 0;
		rows--;
	}
	if (m_rows != 0 && m_rows != rows)
	{
		delete [] p;
		m_rows = 0;
	}
	if (0 == m_rows)
		alloc(m_rows = rows, 1);

	T *out = p;
    int i;
	for(i = 0; i < i0; i++)
		*out++ = v.p[i];
    // Skip row i0
	for(i = i0+1; i < v.m_rows; i++)
		*out++ = v.p[i];
}

/// Paste vector \a v from row i0.
template <typename T>
void vector<T>::paste(int i0, const vector<T>& v)
{
    matrix<T>::paste(i0, 0, v); 
}

} // namespace libNumerics

/// Scalar product.
template <typename T>
T dot(const libNumerics::vector<T>& u, const libNumerics::vector<T>& v)
{
    assert(u.nrow() == v.nrow());
    T d = (T)0;
    for(int i = u.nrow()-1; i >= 0; i--)
        d += u(i)*v(i);
    return d;
}

/// Cross product.
template <typename T>
libNumerics::vector<T> cross(const libNumerics::vector<T>& u,
                             const libNumerics::vector<T>& v)
{
    assert(u.nrow() == 3 && v.nrow() == 3);
    libNumerics::vector<T> w(3);
    w(0) = u(1)*v(2) - u(2)*v(1);
    w(1) = u(2)*v(0) - u(0)*v(2);
    w(2) = u(0)*v(1) - u(1)*v(0);
    return w;
}

/// Sum of all the vector elements
template <typename T>
T sum(const libNumerics::vector<T>& v)
{
	T s = (T)0;
	for (int i = v.nrow()-1; i >= 0; i--)
		s += v(i);
	return s;
}

/// Mean value of the vector
template <typename T>
T mean(const libNumerics::vector<T>& v)
{
	return (T) ( sum(v) / v.nrow() );
}

/// Index of Max value of the vector
template <typename T>
int max(const libNumerics::vector<T>& v)
{
	T m = v[0];
	int idm = 0;
	for (int i = 1; i < v.nrow()-1; i++)
	{
		if (v[i] > m)
		{
			m = v[i];
			idm = i;
		}
	}
	return idm;
}

#endif // ARRAY_H
