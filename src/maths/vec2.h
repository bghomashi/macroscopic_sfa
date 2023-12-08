#pragma once

#include "maths.h"

template <typename T>
struct TVec2 {
    T x, y;

    TVec2& operator+= (const TVec2& o) {
        x += o.x;
        y += o.y;
        return *this;
    }
    TVec2& operator-= (const TVec2& o) {
        x -= o.x;
        y -= o.y;
        return *this;
    }
    TVec2& operator*= (const TVec2& o) {
        x *= o.x;
        y *= o.y;
        return *this;
    }
    TVec2& operator/= (const TVec2& o) {
        x /= o.x;
        y /= o.y;
        return *this;
    }


    TVec2 operator+ (const TVec2& o) const {
        TVec2 r = *this;
        r += o;
        return r;
    }
    TVec2 operator- (const TVec2& o) const {
        TVec2 r = *this;
        r -= o;
        return r;
    }
    TVec2 operator* (const TVec2& o) const {
        TVec2 r = *this;
        r *= o;
        return r;
    }
    TVec2 operator/ (const TVec2& o) const {
        TVec2 r = *this;
        r /= o;
        return r;
    }


    TVec2& operator+= (const T& o) {
        x += o;
        y += o;
        return *this;
    }
    TVec2& operator-= (const T& o) {
        x -= o;
        y -= o;
        return *this;
    }
    TVec2& operator*= (const T& o) {
        x *= o;
        y *= o;
        return *this;
    }
    TVec2& operator/= (const T& o) {
        x /= o;
        y /= o;
        return *this;
    }
    TVec2 operator+ (const T& o) const {
        TVec2 r = *this;
        r += o;
        return r;
    }
    TVec2 operator- (const T& o) const {
        TVec2 r = *this;
        r -= o;
        return r;
    }
    TVec2 operator* (const T& o) const {
        TVec2 r = *this;
        r *= o;
        return r;
    }
    TVec2 operator/ (const T& o) const {
        TVec2 r = *this;
        r /= o;
        return r;
    }

    TVec2& operator=(const T& o) {
        x = o;
        y = o;
        return *this;
    }
};



typedef TVec2<double> dvec2;
typedef TVec2<complex> cvec2;
typedef std::vector<dvec2> d2vector;
typedef std::vector<cvec2> c2vector;

inline double Dot(const dvec2& a, const dvec2& b) {
    return a.x * b.x + a.y * b.y;
};
inline complex Dot(const dvec2& a, const cvec2& b) {
    return a.x * std::conj(b.x) + a.y * std::conj(b.y);
};
inline complex Dot(const cvec2& a, const dvec2& b) {
    return a.x * b.x + a.y * b.y;
};
inline complex Dot(const cvec2& a, const cvec2& b) {
    return a.x * std::conj(b.x) + a.y * std::conj(b.y);
};
inline cvec2 conj(const cvec2& o) {
    return { std::conj(o.x), std::conj(o.y) };
};
inline dvec2 real(const cvec2& o) {
    return { std::real(o.x), std::real(o.y) };
};
inline dvec2 imag(const cvec2& o) {
    return { std::imag(o.x), std::imag(o.y) };
};


inline dvec2 operator+ (const double& o, const dvec2& a) {
    return a + o;
}
inline dvec2 operator* (const double& o, const dvec2& a) {
    return a * o;
}
inline cvec2 operator+ (const double& o, const cvec2& a) {
    return a + o;
}
inline cvec2 operator* (const double& o, const cvec2& a) {
    return a * o;
}
inline cvec2 operator+ (const complex& o, const cvec2& a) {
    return a + o;
}
inline cvec2 operator* (const complex& o, const cvec2& a) {
    return a * o;
}
inline cvec2 operator+ (const complex& o, const dvec2& a) {
    cvec2 r = { a.x, a.y };
    r += o;
    return r;
}
inline cvec2 operator* (const complex& o, const dvec2& a) {
    cvec2 r = { a.x, a.y };
    r *= o;
    return r;
}

inline double length(const dvec2& a) {
    return sqrt(Dot(a, a));
}
inline double length(const cvec2& a) {
    return sqrt(std::real(Dot(a, a)));
}
inline void normalize(dvec2& a) {
    double l = sqrt(Dot(a, a));
    a /= l;
}
inline void normalize(cvec2& a) {
    double l = sqrt(std::real(Dot(a, a)));
    a /= l;
}