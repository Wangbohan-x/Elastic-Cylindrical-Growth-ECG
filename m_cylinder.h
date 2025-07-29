#pragma once
#include"utils.h"

class  Vec3f {

public:
    typedef float ScalarType;
    enum { Dim = 3 };

    Vec3f() {}
    Vec3f(float x, float y, float z) { setValue(x, y, z); }


    explicit Vec3f(const float v[3])
    {
        setValue(v);
    }

    Vec3f& setValue(const float v[3])
    {
        for (int i = 0; i < 3; i++)
            vec[i] = v[i];
        return *this;
    }

    const float* getValue() const
    {
        return vec;
    }

    Vec3f& setValue(float x, float y, float z)
    {
        vec[0] = x;  vec[1] = y;  vec[2] = z;
        return *this;
    }

    void getValue(float& x, float& y, float& z) const
    {
        x = vec[0];  y = vec[1];  z = vec[2];
    }


    operator const float* () const
    {
        return vec;
    }

    operator float* ()
    {
        return vec;
    }

    float X() const { return vec[0]; }
    float Y() const { return vec[1]; }
    float Z() const { return vec[2]; }

    float& X() { return vec[0]; }
    float& Y() { return vec[1]; }
    float& Z() { return vec[2]; }

    /*float& operator [] (int i)
    {
       return vec[i];
    }

    const float& operator [] (int i) const
    {
       return vec[i];
    }*/

    float dot(const Vec3f& v) const
    {
        float s = vec[0] * v.vec[0];
        for (int i = 1; i < 3; i++)
            s += vec[i] * v.vec[i];
        return s;
    }

    float dot(const float* v) const
    {
        return vec[0] * v[0] + vec[1] * v[1] + vec[2] * v[2];
    }

    Vec3f cross(const Vec3f& v) const
    {
        return Vec3f(vec[1] * v.vec[2] - vec[2] * v.vec[1],
            vec[2] * v.vec[0] - vec[0] * v.vec[2],
            vec[0] * v.vec[1] - vec[1] * v.vec[0]);
    }


    float sqrLength() const
    {
        return dot(*this);
    }

    float length() const
    {
        return (float)std::sqrt(dot(*this));
    }

    float normalize()
    {
        float len = length();
        if (len > 0)
            *this /= len;
        return len;
    }

    bool equals(const Vec3f& v, float tolerance) const
    {
        return ((*this - v).sqrLength() <= tolerance);
    }

    friend bool operator == (const Vec3f& v1, const Vec3f& v2)
    {
        for (int i = 0; i < 3; i++)
            if (v1.vec[i] != v2.vec[i])
                return false;
        return true;
    }

    friend bool operator < (const Vec3f& v1, const Vec3f& v2)
    {
        for (int i = 0; i < 3; i++)
            if (v1.vec[i] >= v2.vec[i])
                return false;
        return true;
    }
    friend bool operator <= (const Vec3f& v1, const Vec3f& v2)
    {
        for (int i = 0; i < 3; i++)
            if (v1.vec[i] > v2.vec[i])
                return false;
        return true;
    }
    friend bool operator > (const Vec3f& v1, const Vec3f& v2)
    {
        for (int i = 0; i < 3; i++)
            if (v1.vec[i] <= v2.vec[i])
                return false;
        return true;
    }
    friend bool operator >= (const Vec3f& v1, const Vec3f& v2)
    {
        for (int i = 0; i < 3; i++)
            if (v1.vec[i] < v2.vec[i])
                return false;
        return true;
    }

    friend bool operator != (const Vec3f& v1, const Vec3f& v2)
    {
        return !(v1 == v2);
    }

    Vec3f& operator *= (float s)
    {
        for (int i = 0; i < 3; i++)
            vec[i] *= s;
        return *this;
    }

    Vec3f& operator /= (float s)
    {
        for (int i = 0; i < 3; i++)
            vec[i] /= s;
        return *this;
    }

    Vec3f& operator += (const Vec3f& v)
    {
        for (int i = 0; i < 3; i++)
            vec[i] += v.vec[i];
        return *this;
    }

    Vec3f& operator -= (const Vec3f& v)
    {
        for (int i = 0; i < 3; i++)
            vec[i] -= v.vec[i];
        return *this;
    }

    void negate()
    {
        for (int i = 0; i < 3; i++)
            vec[i] = -vec[i];
    }


    friend Vec3f operator * (float s, const Vec3f& v) { return Vec3f(v) *= s; }
    friend Vec3f operator * (const Vec3f& v, float s) { return Vec3f(v) *= s; }
    friend Vec3f operator / (const Vec3f& v, float s) { return Vec3f(v) /= s; }
    friend Vec3f operator + (const Vec3f& v1, const Vec3f& v2) { return Vec3f(v1) += v2; }
    friend Vec3f operator - (const Vec3f& v1, const Vec3f& v2) { return Vec3f(v1) -= v2; }
    friend Vec3f operator - (const Vec3f& v1) { Vec3f v2 = v1; v2.negate(); return v2; }
    Vec3f operator*(const Vec3f& v) const { return Vec3f(vec[0] * v.vec[0], vec[1] * v.vec[1], vec[2] * v.vec[2]); }

    // Konvertierung:

    Vec3f& setValue(const Vec3f& v)
    {
        for (int i = 0; i < 3; i++)
            vec[i] = (float)v[i];
        return *this;
    }

    void getValue(Vec3f& v) const
    {
        for (int i = 0; i < 3; i++)
            v[i] = (float)vec[i];
    }

private:
    float vec[3];
};

struct CylinderElement
{
	double radius;
    std::vector<Vec3f> positions;//Ô²ÖùÌå¶Ëµã
    Vec3f direction;
};

double GetHeight(Point_set& points, double dx, double dy, double dz, double cx, double cy, double cz) {
	double min_proj = DBL_MAX;
	double max_proj = DBL_MIN;
	double proj = 0;
	for (size_t i = 0; i < points.size(); i++)
	{
		Point point = points.point(i);
		Point vec(point.x() - cx, point.y() - cy, point.z() - cz);
		proj = vec.x() * dx + vec.y() * dy + vec.z() * dz;
		if (proj < min_proj) min_proj = proj;
		if (proj > max_proj) max_proj = proj;
	}
	return max_proj - min_proj;
}
