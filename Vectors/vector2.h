#pragma once
#include <cmath>
#include <iostream>

template <typename T>
class Vector2 {
public:
    T x, y;

    Vector2() : x(0), y(0) {}
    Vector2(T x, T y) : x(x), y(y) {}

    Vector2 operator+(const Vector2& other) const {
        return Vector2(x + other.x, y + other.y);
    }

    Vector2& operator+=(const Vector2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2 operator-(const Vector2& other) const {
        return Vector2(x - other.x, y - other.y);
    }

    Vector2& operator-=(const Vector2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2 operator*(T scalar) const {
        return Vector2(x * scalar, y * scalar);
    }

    Vector2& operator*=(T scalar) {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2 operator/(T scalar) const {
        return Vector2(x / scalar, y / scalar);
    }

    Vector2& operator/=(T scalar) {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    T dot(const Vector2& other) const {
        return x * other.x + y * other.y;
    }

    T sqlength() const {
        return x * x + y * y;
    }

    T length() const {
        return std::sqrt(sqlength());
    }

    Vector2 normalized() const {
        T len = length();
        if (len == 0) return Vector2(0, 0);
        return *this / len;
    }

    void normalize() {
        T len = length();
        if (len != 0) {
            *this /= len;
        }
    }

    friend Vector2 operator*(T scalar, const Vector2& vec) {
        return Vector2(vec.x * scalar, vec.y * scalar);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector2& vec) {
        os << "(" << vec.x << ", " << vec.y << ")";
        return os;
    }
};