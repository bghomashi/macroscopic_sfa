#pragma once

#include "maths.h"

template <typename T>
struct TVec3 {
    T x, y, z;

    TVec3& operator+= (const TVec3& o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    TVec3& operator-= (const TVec3& o) {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }
    TVec3& operator*= (const TVec3& o) {
        x *= o.x;
        y *= o.y;
        z *= o.z;
        return *this;
    }
    TVec3& operator/= (const TVec3& o) {
        x /= o.x;
        y /= o.y;
        z /= o.z;
        return *this;
    }


    TVec3 operator+ (const TVec3& o) const {
        TVec3 r = *this;
        r += o;
        return r;
    }
    TVec3 operator- (const TVec3& o) const {
        TVec3 r = *this;
        r -= o;
        return r;
    }
    TVec3 operator* (const TVec3& o) const {
        TVec3 r = *this;
        r *= o;
        return r;
    }
    TVec3 operator/ (const TVec3& o) const {
        TVec3 r = *this;
        r /= o;
        return r;
    }


    TVec3& operator+= (const T& o) {
        x += o;
        y += o;
        z += o;
        return *this;
    }
    TVec3& operator-= (const T& o) {
        x -= o;
        y -= o;
        z -= o;
        return *this;
    }
    TVec3& operator*= (const T& o) {
        x *= o;
        y *= o;
        z *= o;
        return *this;
    }
    TVec3& operator/= (const T& o) {
        x /= o;
        y /= o;
        z /= o;
        return *this;
    }
    TVec3 operator+ (const T& o) const {
        TVec3 r = *this;
        r += o;
        return r;
    }
    TVec3 operator- (const T& o) const {
        TVec3 r = *this;
        r -= o;
        return r;
    }
    TVec3 operator* (const T& o) const {
        TVec3 r = *this;
        r *= o;
        return r;
    }
    TVec3 operator/ (const T& o) const {
        TVec3 r = *this;
        r /= o;
        return r;
    }

    TVec3& operator=(const T& o) {
        x = o;
        y = o;
        return *this;
    }
};



typedef TVec3<double> dvec3;
typedef TVec3<complex> cvec3;
typedef std::vector<dvec3> d3vector;
typedef std::vector<cvec3> c3vector;

inline double Dot(const dvec3& a, const dvec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
};
inline complex Dot(const dvec3& a, const cvec3& b) {
    return a.x * std::conj(b.x) + a.y * std::conj(b.y) + a.z * std::conj(b.z);
};
inline complex Dot(const cvec3& a, const dvec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
};
inline complex Dot(const cvec3& a, const cvec3& b) {
    return a.x * std::conj(b.x) + a.y * std::conj(b.y) + a.z * std::conj(b.z);
};
inline cvec3 conj(const cvec3& o) {
    return { std::conj(o.x), std::conj(o.y), std::conj(o.z) };
};
inline dvec3 real(const cvec3& o) {
    return { std::real(o.x), std::real(o.y), std::real(o.z) };
};
inline dvec3 imag(const cvec3& o) {
    return { std::imag(o.x), std::imag(o.y), std::imag(o.z) };
};


inline dvec3 operator+ (const double& o, const dvec3& a) {
    return a + o;
}
inline dvec3 operator* (const double& o, const dvec3& a) {
    return a * o;
}
inline cvec3 operator+ (const double& o, const cvec3& a) {
    return a + o;
}
inline cvec3 operator* (const double& o, const cvec3& a) {
    return a * o;
}
inline cvec3 operator+ (const complex& o, const cvec3& a) {
    return a + o;
}
inline cvec3 operator* (const complex& o, const cvec3& a) {
    return a * o;
}
inline cvec3 operator+ (const complex& o, const dvec3& a) {
    cvec3 r = { a.x, a.y, a.z };
    r += o;
    return r;
}
inline cvec3 operator* (const complex& o, const dvec3& a) {
    cvec3 r = { a.x, a.y, a.z };
    r *= o;
    return r;
}

inline double length(const dvec3& a) {
    return sqrt(Dot(a, a));
}
inline double length(const cvec3& a) {
    return sqrt(std::real(Dot(a, a)));
}
inline void normalize(dvec3& a) {
    double l = sqrt(Dot(a, a));
    a /= l;
}
inline void normalize(cvec3& a) {
    double l = sqrt(std::real(Dot(a, a)));
    a /= l;
}