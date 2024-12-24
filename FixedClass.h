//
// Created by rcmpl on 23.12.2024.
//

#ifndef CPP_PROJECT_2_FIXEDCLASS_H
#define CPP_PROJECT_2_FIXEDCLASS_H
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <compare>
namespace types {
    template<size_t N, size_t K, bool Fast = false>
    class Fixed {
        using ValueType = typename std::conditional_t<
                Fast,
                typename std::conditional_t<
                        (N <= 8), int_fast8_t,
                        typename std::conditional_t<
                                (N <= 16), int_fast16_t,
                                typename std::conditional_t<
                                        (N <= 32), int_fast32_t,
                                        typename std::conditional_t<
                                                (N <= 64), int_fast64_t,
                                                void
                                        >
                                >
                        >
                >,
                typename std::conditional_t<
                        (N == 8), int8_t,
                        typename std::conditional_t<
                                (N == 16), int16_t,
                                typename std::conditional_t<
                                        (N == 32), int32_t,
                                        typename std::conditional_t<
                                                (N == 64), int64_t,
                                                void
                                        >
                                >
                        >
                >
        >;
    public:
        static constexpr std::size_t NFromFixed = N;
        static constexpr std::size_t KFromFixed = K;
        static constexpr bool FastFromFixed = Fast;

        constexpr Fixed() : value(0) {}

        constexpr explicit Fixed(int32_t value) : value(value << K) {}

        constexpr explicit Fixed(double value) : value(value * (1 << K)) {}

        constexpr explicit Fixed(float value) : value(value * (1 << K)) {}

        static constexpr Fixed from_raw(ValueType x) {
            Fixed ret;
            ret.value = x;
            return ret;
        }

        explicit operator double() const noexcept {
            return static_cast<double>(this->value) / (1 << K);
        }

        explicit operator float() const noexcept {
            return static_cast<float>(this->value) / (1 << K);
        }

        auto operator<=>(const Fixed &) const = default;

        bool operator==(const Fixed &) const = default;

        template<size_t N1, size_t K1, bool Fast1>
        friend std::ostream &operator<<(std::ostream &os, const Fixed<N1, K1, Fast1> &f);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator+(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator-(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator*(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator/(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator/(const double& a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator*(const Fixed<N1, K1, Fast1> &a, const double &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator*(const double& a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator <(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator >(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator <=(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator >=(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        template<size_t N1, size_t K1, bool Fast1>
        friend Fixed<N1, K1, Fast1> operator ==(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b);

        ValueType value;
    };

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator==(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value == b.value;
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator>=(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value >= b.value;
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator<=(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value <= b.value;
    }
    template <size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> min(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value < b.value ? a : b;
    }
    double min(double a, double b) {
        return a < b ? a : b;
    }
    float min(float a, float b) {
        return a < b ? a : b;
    }
    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator>(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value > b.value;
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator<(const Fixed<N1, K1, Fast1> &a, const Fixed<N1, K1, Fast1> &b) {
        return a.value < b.value;
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator*(const double &a, const Fixed<N1, K1, Fast1> &b) {
        return Fixed<N1, K1, Fast1>::from_raw((a * b.value) >> K1);
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator*(const Fixed<N1, K1, Fast1> &a, const double &b) {
        return Fixed<N1, K1, Fast1>::from_raw((a.value * b)  >> K1);
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator/(const double &a, const Fixed<N1, K1, Fast1> &b) {
        return Fixed<N1, K1, Fast1>::from_raw((a * (1 << K1)) / b.value);
    }

    template<size_t N1, size_t K1, bool Fast1>
    Fixed<N1, K1, Fast1> operator/(const Fixed<N1, K1, Fast1> &a, const double &b) {
        return Fixed<N1, K1, Fast1>::from_raw((a.value << K1) / b);
    }

    template<size_t N, size_t K, bool Fast>
    std::istream &operator>>(std::istream &is, Fixed<N, K, Fast> &f) {
        int x;
        is >> x;
        f = x;
        return is;
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> operator+(const Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return Fixed<N, K, Fast>::from_raw(a.value + b.value);
    }

    template<std::size_t N, std::size_t K>
    using FastFixed = Fixed<N, K, true>;

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> operator-(const Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return Fixed<N, K, Fast>::from_raw(a.value - b.value);
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> operator*(const Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return Fixed<N, K, Fast>::from_raw((int64_t(a.value) * b.value) >> K);
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> operator/(const Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return Fixed<N, K, Fast>::from_raw((a.value << K) / b.value);
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> &operator+=(Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return a = a + b;
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> &operator-=(Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return a = a - b;
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> &operator*=(Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return a = a * b;
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> &operator/=(Fixed<N, K, Fast> &a, const Fixed<N, K, Fast> &b) {
        return a = a / b;
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> operator-(const Fixed<N, K, Fast> &a) {
        return Fixed<N, K, Fast>::from_raw(-a.value);
    }

    template<size_t N, size_t K, bool Fast>
    std::ostream &operator<<(std::ostream &os, const Fixed<N, K, Fast> &f) {
        return os << static_cast<double>(f.value) / (1 << K);
    }

    template<size_t N, size_t K, bool Fast>
    Fixed<N, K, Fast> abs(const Fixed<N, K, Fast> &a) {
        return a.value < 0 ? -a : a;
    }
}
#endif //CPP_PROJECT_2_FIXEDCLASS_H
