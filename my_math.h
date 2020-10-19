#pragma once


#include <cmath>
#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>


template<int64_t Columns, int64_t Rows, typename T>
class Matrix;

template <int64_t Dimension, typename T>
struct Vector
{
	Vector()
	{
		for (int64_t i = 0; i < Dimension; i++)
			s_data[i] = T();
	}

	T &operator[] (const int64_t i)
	{
		assert(i < Dimension);
		return s_data[i];
	}

	const T &operator[] (const int64_t i) const
	{
		assert(i < Dimension);
		return s_data[i];
	}

private:
	T s_data[Dimension];
};

template <typename T>
struct Vector<2, T>
{
	Vector()
	: s_x(T()), s_y(T())
	{}

	Vector(T x, T y)
	: s_x(x), s_y(y)
	{}

	template <typename U>
	Vector<2, T> (const Vector<2, U> &vector);

	T &operator[] (const int64_t i)
	{
		assert(i < 2);
		return i <= 0 ? s_x : s_y;
	}

	const T &operator[] (const int64_t i) const
	{
		assert(i < 2);
		return i <= 0 ? s_x : s_y;
	}

	T s_x, s_y;
};

template <typename T>
struct Vector<3, T>
{
	Vector()
	: s_x(T()), s_y(T()), s_z(T())
	{ }

	Vector(T x, T y, T z)
	: s_x(x), s_y(y), s_z(z)
	{ }

	template <typename U>
	Vector<3, T> (const Vector<3, U> &vector);

	T &operator[] (const int64_t i)
	{
		assert(i < 3);
		return i <= 0 ? s_x : (1 == i ? s_y : s_z);
	}

	const T &operator[](const int64_t i) const
	{
		assert(i < 3);
		return i <= 0 ? s_x : (1 == i ? s_y : s_z);
	}

	float norm()
	{
		return std::sqrt(s_x * s_x + s_y * s_y + s_z * s_z);
	}

	Vector<3, T> &normalize(T length = 1)
	{
		length = length * (length / norm());
		return length;
	}

	T s_x, s_y, s_z;
};

template<int64_t Dimension, typename T>
T operator* (const Vector<Dimension, T> &left_hand_side, const Vector<Dimension, T> &right_hand_side)
{
	T result = T();

	for (int64_t i = 0; i < Dimension; i++)
		result += left_hand_side[i] * right_hand_side[i];

	return result;
}

template<int64_t Dimension, typename T>
Vector<Dimension, T> operator+ (Vector<Dimension, T> left_hand_side, const Vector<Dimension, T> &right_hand_side)
{
	for (uint64_t i = 0; i < Dimension; i++)
		left_hand_side[i] += right_hand_side[i];

	return left_hand_side;
}

template<int64_t Dimension, typename T>
Vector<Dimension, T> operator- (Vector<Dimension, T> left_hand_side, const Vector<Dimension, T> &right_hand_side)
{
	for (int64_t i = 0; i < Dimension; i++)
		left_hand_side[i] -= right_hand_side[i];

	return left_hand_side;
}

template<int64_t Dimension, typename T, typename U>
Vector<Dimension, T> operator* (Vector<Dimension, T> left_hand_side, const U &right_hand_side)
{
	for (int64_t i = 0; i < Dimension; i++)
		left_hand_side[i] *= right_hand_side;

	return left_hand_side;
}

template<int64_t Dimension, typename T, typename U>
Vector<Dimension, T> operator/ (Vector<Dimension, T> left_hand_side, const U &right_hand_side)
{
	for (int64_t i = 0; i < Dimension; i++)
		left_hand_side[i] /= right_hand_side;

	return left_hand_side;
}

template<int64_t Length, int64_t Dimension, typename T>
Vector<Length, T> embed(const Vector<Dimension, T> &vector, T fill = 1)
{
	Vector<Length, T> result;

	for (int64_t i = 0; i < Length; i++)
		result[i] = (i < Dimension ? vector[i] : fill);

	return result;
}

template<int64_t Length, int64_t Dimension, typename T>
Vector<Length, T> projection(const Vector<Dimension, T> &vector)
{
	Vector<Length, T> result;

	for (int64_t i = 0; i < Length; i++)
		result[i] = vector[i];

	return result;
}

template <typename T>
Vector<3, T> cross(Vector<3, T> first_vector, Vector<3, T> second_vector)
{
	return Vector<3, T>(first_vector.s_y * second_vector.s_z - first_vector.s_z * second_vector.s_y,
					    first_vector.s_z * second_vector.s_x - first_vector.s_x * second_vector.s_z,
						first_vector.s_x * second_vector.s_y - first_vector.s_y * second_vector.s_x);
}

template <int64_t Dimension, typename T>
std::ostream &operator<< (std::ostream &output, Vector<Dimension, T> &vector)
{
	for(int64_t i = 0; i < Dimension; i++)
		output << vector[i] << " " ;

	return output ;
}

template<int64_t Dimension, typename T>
struct dt
{
	static T det(const Matrix<Dimension, Dimension, T> &source)
	{
		T result = 0;

		for (int64_t i = 0; i < Dimension; i ++)
			result += source[0][i] * source.cofactor(0, i);

		return result;
	}
};

template<typename T>
struct dt<1, T>
{
	static T det(const Matrix<1, 1, T> &source)
	{
		return source[0][0];
	}
};

template<int64_t Rows, int64_t Columns, typename T>
class Matrix
{
public:
	Matrix() = default;

	Vector<Columns, T> &operator[] (const int64_t i)
	{
		assert(i < Rows);
		return rows[i];
	}

	const Vector<Columns, T> &operator[] (const int64_t i) const
	{
		assert(i < Rows);
		return rows[i];
	}

	Vector<Rows, T> column(const int64_t j) const
	{
		assert(j < Columns);
		Vector<Rows, T> result;

		for (int64_t i = 0; i < Rows; i++)
			result[i] = rows[i][j];

		return result;
	}

	void set_column(int64_t j, Vector < Rows, T> vector)
	{
		assert(j < Columns);

		for (int64_t i = 0; i < Rows; i++)
			rows[i][j] = vector[i];
	}

	static Matrix<Rows, Columns, T>
	identity()
	{
		Matrix<Rows, Columns, T> result;

		for (int64_t i = 0; i < Rows; i++)
			for (int64_t j = 0; j < Columns; j++)
				result[i][j] = (i == j);

		return result;
	}

	T det() const
	{
		return dt<Columns, T>::det(*this);
	}

	Matrix<Rows - 1, Columns - 1, T>
	get_minor(int64_t row, int64_t column) const
	{
		Matrix<Rows - 1, Columns - 1, T> result;

		for (int64_t i = 0; i < Rows - 1; i++)
			for (int64_t j = 0; j < Columns - 1; j++)
				result[i][j] = rows[i < row ? i : i + 1][j < column ? j : j + 1];

		return result;
	}

	T cofactor(int64_t row, int64_t column) const
	{
		return get_minor(row, column).det() * ((row + column) % 2 ? -1 : 1);
	}

	Matrix<Rows, Columns, T>
	adjugate_matrix() const
	{
		Matrix<Rows ,Columns, T> result;

		for (int64_t i = 0; i < Rows; i++)
			for (int64_t j = 0; j < Columns; j++)
				result[i][j] = cofactor(i, j);

		return result;
	}

	Matrix<Rows, Columns, T>
	invert_transpose()
	{
		Matrix<Rows ,Columns, T> result = adjugate_matrix();
		T temp = result[0] * rows[0];
		return result / temp;
	}

	Matrix<Columns, Rows, T>
	invert()
	{
		return invert_transpose().transpose();
	}

	Matrix<Columns, Rows, T>
	transpose()
	{
		Matrix<Columns, Rows, T> result;

		for (int64_t i = 0; i < Rows; i++)
			result[i] = column(i);

		return result;
	}

private:
	Vector<Columns, T> rows[Rows];
};

template<size_t Rows, size_t Columns, typename T>
Vector<Rows, T> operator* (const Matrix<Rows, Columns, T> &left_hand_side, const Vector<Columns, T> &right_hand_side)
{
	Vector<Rows, T> result;

	for (int64_t i = 0; i < Rows; i++)
		result[i] = left_hand_side[i] * right_hand_side;

	return result;
}

template<int64_t R1, int64_t C1, int64_t C2, typename T>
Matrix<R1, C2, T> operator* (const Matrix<R1, C1, T> &left_hand_side, const Matrix<C1, C2, T> &right_hand_side)
{
	Matrix<R1, C2, T> result;

	for (int64_t i = 0; i < R1; i++)
		for (int64_t j = 0; j < C2; j++)
			result[i][j] = left_hand_side[i] * right_hand_side.column(j);

	return result;
}

template<int64_t Rows, int64_t Columns, typename T>
Matrix<Columns, Rows, T> operator/ (Matrix<Rows, Columns, T> left_hand_side, const T &right_hand_side)
{
	for (int64_t i = 0; i < Rows; i++)
		left_hand_side[i] = left_hand_side[i] / right_hand_side;

	return left_hand_side;
}

template<int64_t Rows, int64_t Columns, typename T>
std::ostream &operator<< (std::ostream &output, Matrix<Rows, Columns, T> &matrix)
{
	for (int64_t i = 0; i < Rows; i++)
		output << matrix[i] << std::endl;

	return output;
}

using VectorTwoInt     = Vector<2,    int>;
using VectorThreeInt   = Vector<3,    int>;
using VectorTwoFloat   = Vector<2,    float>;
using VectorThreeFloat = Vector<3,    float>;
using VectorFourFloat  = Vector<4,    float>;
using MatrixFourByFour = Matrix<4, 4, float>;