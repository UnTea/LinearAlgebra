#pragma once


#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>


template<size_t dimension_columns, size_t dimension_rows, typename Type> class Matrix;

template <size_t dimension, typename Type>
struct Vector
{
	Vector()
	{ for (size_t i(dimension); i--; s_data[i] = Type()); }

	Type &operator[](const size_t index)
	{
		assert(index < dimension);
		return s_data[index];
	}

	const Type &operator[](const size_t index) const
	{
		assert(index < dimension);
		return s_data[index];
	}

private:
	Type s_data[dimension];
};

template <typename Type>
struct Vector<2, Type>
{
	Vector()
	: s_x(Type()), s_y(Type())
	{ }

	Vector(Type x, Type y)
	: s_x(x), s_y(y)
	{ }

	template <typename U> explicit Vector<2, Type> (const Vector<2, U> &vector);

	Type &operator[](const size_t index)
	{
		assert(index < 2);
		return index <= 0 ? s_x : s_y;
	}

	const Type &operator[](const size_t index) const
	{
		assert(index < 2);
		return index <= 0 ? s_x : s_y;
	}

	Type s_x, s_y;
};

template <typename Type>
struct Vector<3, Type>
{
	Vector()
	: s_x(Type()), s_y(Type()), s_z(Type())
	{ }

	Vector(Type x, Type y, Type z)
	: s_x(x), s_y(y), s_z(z)
	{ }

	template <typename U> explicit Vector<3, Type>(const Vector<3, U> &vector);

	Type &operator[](const size_t index)
	{
		assert(index < 3);
		return index <= 0 ? s_x : (1 == index ? s_y : s_z);
	}

	const Type &operator[](const size_t index) const
	{
		assert(index < 3);
		return index <= 0 ? s_x : (1 == index ? s_y : s_z);
	}

	float norm()
	{
		return std::sqrt(s_x * s_x + s_y * s_y + s_z * s_z);
	}

	Vector<3, Type> &normalize(Type length = 1)
	{
		length = length * (length / norm());
		return length;
	}

	Type s_x, s_y, s_z;
};

template<size_t dimension, typename Type>
Type operator*(const Vector<dimension, Type> &left_hand_side, const Vector<dimension, Type> &right_hand_side)
{
	Type result = Type();
	for (size_t i(dimension); i--; result += left_hand_side[i] * right_hand_side[i]);
	return result;
}

template<size_t dimension, typename Type>
Vector<dimension, Type> operator+(Vector<dimension, Type> left_hand_side, const Vector<dimension, Type> &right_hand_side)
{
	for (size_t i(dimension); i--; left_hand_side[i] += right_hand_side[i]);
	return left_hand_side;
}

template<size_t dimension, typename Type>
Vector<dimension, Type> operator-(Vector<dimension, Type> left_hand_side, const Vector<dimension, Type> &right_hand_side)
{
	for (size_t i(dimension); i--; left_hand_side[i] -= right_hand_side[i]);
	return left_hand_side;
}

template<size_t Dimension, typename Type, typename U>
Vector<Dimension, Type> operator*(Vector<Dimension,Type> left_hand_side, const U &right_hand_side)
{
	for (size_t i(Dimension); i--; left_hand_side[i] *= right_hand_side);
	return left_hand_side;
}

template<size_t Dimension, typename Type, typename U>
Vector<Dimension, Type> operator/(Vector<Dimension, Type> left_hand_side, const U &right_hand_side)
{
	for (size_t i(Dimension); i--; left_hand_side[i] /= right_hand_side);
	return left_hand_side;
}

template<size_t length, size_t dimension, typename Type>
Vector<length, Type> embed(const Vector<dimension, Type> &vector, Type fill = 1)
{
	Vector<length, Type> result;
	for (size_t i(length); i--; result[i] = (i < dimension ? vector[i] : fill));
	return result;
}

template<size_t length, size_t dimension, typename Type>
Vector<length, Type> projection(const Vector<dimension, Type> &vector)
{
	Vector<length, Type> result;
	for (size_t i(length); i--; result[i] = vector[i]);
	return result;
}

template <typename Type> Vector<3, Type> cross(Vector<3, Type> first_vector, Vector<3, Type> second_vector)
{
	return Vector<3, Type>(first_vector.s_y * second_vector.s_z - first_vector.s_z * second_vector.s_y,
						   first_vector.s_z * second_vector.s_x - first_vector.s_x * second_vector.s_z,
						   first_vector.s_x * second_vector.s_y - first_vector.s_y * second_vector.s_x);
}

template <size_t dimension, typename Type>
std::ostream &operator<<(std::ostream &output, Vector<dimension, Type> &vector)
{
	for(unsigned int i(0); i < dimension; i++)
		output << vector[i] << " " ;
	return output ;
}

template<size_t dimension, typename Type>
struct dt
{
	static Type det(const Matrix<dimension, dimension, Type> &source)
	{
		Type result(0);
		for (size_t i(dimension); i--; result += source[0][i] * source.cofactor(0, i));
		return result;
	}
};

template<typename Type>
struct dt<1, Type>
{
	static Type det(const Matrix<1, 1, Type> &source)
	{ return source[0][0]; }
};

template<size_t dimension_rows, size_t dimension_columns, typename Type>
class Matrix
{
	Vector<dimension_columns, Type> rows[dimension_rows];

public:
	Matrix() = default;

	Vector<dimension_columns, Type> &operator[](const size_t index)
	{
		assert(index < dimension_rows);
		return rows[index];
	}

	const Vector<dimension_columns, Type> &operator[](const size_t index) const
	{
		assert(index < dimension_rows);
		return rows[index];
	}

	Vector<dimension_rows, Type> column(const size_t index) const
	{
		assert(index < dimension_columns);
		Vector<dimension_rows, Type> result;
		for (size_t i(dimension_rows); i--; result[i] = rows[i][index]);
		return result;
	}

	void set_column(size_t j, Vector < dimension_rows, Type> vector)
	{
		assert(j < dimension_columns);
		for (size_t i(dimension_rows); i--; rows[i][j] = vector[i]);
	}

	static Matrix<dimension_rows, dimension_columns, Type> identity()
	{
		Matrix<dimension_rows, dimension_columns, Type> result;
		for (size_t i(dimension_rows); i--; )
			for (size_t j(dimension_columns); j--; result[i][j] = (i == j));
		return result;
	}

	Type det() const
	{
		return dt<dimension_columns, Type>::det(*this);
	}

	Matrix<dimension_rows - 1, dimension_columns - 1, Type> get_minor(size_t row, size_t column) const
	{
		Matrix<dimension_rows - 1, dimension_columns - 1, Type> result;
		for (size_t i(dimension_rows - 1); i--; )
			for (size_t j(dimension_columns - 1); j--; result[i][j] = rows[i < row ? i : i + 1][j < column ? j : j + 1]);
		return result;
	}

	Type cofactor(size_t row, size_t column) const
	{
		return get_minor(row, column).det() * ((row + column) % 2 ? -1 : 1);
	}

	Matrix<dimension_rows, dimension_columns, Type> adjugate_matrix() const
	{
		Matrix<dimension_rows ,dimension_columns, Type> result;
		for (size_t i(dimension_rows); i--; )
			for (size_t j(dimension_columns); j--; result[i][j] = cofactor(i, j));
		return result;
	}

	Matrix<dimension_rows, dimension_columns, Type> invert_transpose()
	{
		Matrix<dimension_rows ,dimension_columns, Type> result = adjugate_matrix();
		Type temp = result[0] * rows[0];
		return result / temp;
	}

	Matrix<dimension_columns, dimension_rows, Type> invert()
	{ return invert_transpose().transpose(); }

	Matrix<dimension_columns, dimension_rows, Type> transpose()
	{
		Matrix<dimension_columns, dimension_rows, Type> result;
		for (size_t i(dimension_rows); i--; result[i] = column(i));
		return result;
	}
};

template<size_t dimension_rows, size_t dimension_columns, typename Type>
Vector<dimension_rows, Type> operator*(const Matrix<dimension_rows, dimension_columns, Type> &left_hand_side, const Vector<dimension_columns, Type> &right_hand_side)
{
	Vector<dimension_rows, Type> result;
	for (size_t i(dimension_rows); i--; result[i] = left_hand_side[i] * right_hand_side);
	return result;
}

template<size_t R1, size_t C1, size_t C2, typename Type>
Matrix<R1, C2, Type> operator*(const Matrix<R1, C1, Type> &left_hand_side, const Matrix<C1, C2, Type> &right_hand_side)
{
	Matrix<R1, C2, Type> result;
	for (size_t i(R1); i--; )
		for (size_t j(C2); j--; result[i][j] = left_hand_side[i] * right_hand_side.column(j));
	return result;
}

template<size_t dimension_rows, size_t dimension_columns, typename Type>
Matrix<dimension_columns, dimension_rows, Type> operator/(Matrix<dimension_rows, dimension_columns, Type> left_hand_side, const Type &right_hand_side)
{
	for (size_t i(dimension_rows); i--; left_hand_side[i] = left_hand_side[i] / right_hand_side);
	return left_hand_side;
}

template<size_t dimension_rows, size_t dimension_columns, class Type>
std::ostream &operator<<(std::ostream &output, Matrix<dimension_rows, dimension_columns, Type> &matrix)
{
	for (size_t i(0); i < dimension_rows; i++)
		output << matrix[i] << std::endl;
	return output;
}

using VectorTwoInt     = Vector<2, int>;
using VectorTwoFloat   = Vector<2, float>;
using VectorThreeInt   = Vector<3, int>;
using VectorThreeFloat = Vector<3, float>;
using VectorFourFloat  = Vector<4, float>;
using MatrixFourByFour = Matrix<4, 4, float>;