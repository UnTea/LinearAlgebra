#include "my_math.h"

template <> template <> Vector<3, int>::Vector(const Vector<3, float> &vector)
: s_x(int(vector.s_x + .5f)), s_y(int(vector.s_y + .5f)), s_z(int(vector.s_z + .5f))
{ }

template <> template <> Vector<3, float>::Vector(const Vector<3, int> &vector)
: s_x(vector.s_x), s_y(vector.s_y), s_z(vector.s_z)
{ }

template <> template <> Vector<2, int>::Vector(const Vector<2, float> &vector)
: s_x(int(vector.s_x + .5f)), s_y(int(vector.s_y + .5f))
{ }

template <> template <> Vector<2, float>::Vector(const Vector<2, int> &vector)
: s_x(vector.s_x), s_y(vector.s_y)
{ }