///
/// (C) 2018 CHIV, All rights reserved.
///
/// \brief 提供一些方便的计算函数.
#ifndef _b_math_h
#define _b_math_h

#include "btypes_t.h"

namespace btypes {
/// bmath提供的函数都需要3维点, 2维点是三维点的特例, 不单独给出.
/// 接受任意类型, 如果是装箱类型, 必须满足给出double出口, 否则不能使用.
/// 因为所有的返回值尽可能向double靠拢. 下面例如.
/// class MyDouble {
/// public:
///      operator double()
///      {
///          return x;
///      }
///      MyDouble(double x_)
///      {
///          x = x_;
///      }
///      void operator=(double x_)
///      {
///          x = x_;
///      }
///      double x;
/// };
namespace bmath {
    template <typename T = double, int M>
    double length(const btypes::Array<T, M>& vec)
    {
        double sum = 0;
        for (int k = 0; k < vec.Size(); ++k) {
            sum += vec[k] * vec[k];
        }
        return sqrt(sum);
    }
    /// 求两点之间的距离
    template <typename T = double, int M>
    double distance(const btypes::Array<T, M>& pt1, const btypes::Array<T, M>& pt2)
    {
        return length(pt2 - pt1);
    }
    /// 求点到直线的距离, 直线用通用式表示 Ax+By+C = 0
    template <typename T = double, int M>
    double distance(const btypes::Array<T, M>& pt, T A, T B, T C)
    {
        return fabs(A * pt[0] + B * pt[1] + C) / sqrt(A * A + B * B);
    }
    /// 求点到面的距离, 面用通用式表示 Ax+By+Cz+D = 0
    template <typename T = double, int M>
    double distance(const btypes::Array<T, M>& pt, T A, T B, T C, T D)
    {
        return fabs(A * pt[0] + B * pt[1] + C * pt[2] + D) / sqrt(A * A + B * B + C * C);
    }
    /// 求点到直线的距离, 带符号
    template <typename T = double, int M>
    double s_distance(const btypes::Array<T, M>& pt, T A, T B, T C)
    {
        return (A * pt[0] + B * pt[1] + C) / sqrt(A * A + B * B);
    }
    /// 求点到面的距离, 带符号
    template <typename T = double, int M>
    double s_distance(const btypes::Array<T, M>& pt, T A, T B, T C, T D)
    {
        return (A * pt[0] + B * pt[1] + C * pt[2] + D) / sqrt(A * A + B * B + C * C);
    }
    /// 求点积
    template <typename T = double, int M>
    double dot(const btypes::Array<T, M>& vec1, const btypes::Array<T, M>& vec2)
    {
        assert(vec1.Size() == vec2.Size());
        btypes::Array<T> v1 = vec1 / length(vec1);
        btypes::Array<T> v2 = vec2 / length(vec2);
        double sum = 0;
        for (int k = 0; k < vec1.Size(); ++k) {
            sum += v1[k] * v2[k];
        }
        return sum;
    }
    // 求叉积, 只有空间向量支持叉乘
    template <typename T = double, int M>
    btypes::Array<T> cross(const btypes::Array<T, M>& vec1, const btypes::Array<T, M>& vec2)
    {
        assert(vec1.Size() == vec2.Size() && vec1.Size() == 3);
        return btypes::Array<T>(
            vec1[1] * vec2[2] - vec1[2] * vec2[1],
            vec1[2] * vec2[0] - vec1[0] * vec2[2],
            vec1[0] * vec2[1] - vec1[1] * vec2[0]);
    }
    // 求点到直线的距离, 直线由其上的两点表示, 如果pt1和pt2为同一点, 将造成除0异常
    template <typename T = double, int M>
    double distance(
        const btypes::Array<T, M>& pt,
        const btypes::Array<T, M>& pt1, const btypes::Array<T, M>& pt2)
    {
        assert(pt.Size() == pt1.Size() && pt.Size() == pt2.Size() && pt.Size() == 3);
        return length(cross(pt - pt1, pt2 - pt1)) / length(pt2 - pt1);
    }
    // 过三点求面方程通用式 Ax+By+Cx+D=0
    // 返回数组 [A,B,C,D]
    template <typename T = double>
    btypes::Array<T> find_equation(
        T x1, T y1, T z1,
        T x2, T y2, T z2,
        T x3, T y3, T z3)
    {
        return btypes::Array<T>(
            y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2,
            -x1 * z2 + x1 * z3 + x2 * z1 - x2 * z3 - x3 * z1 + x3 * z2,
            x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2,
            -x1 * y2 * z3 + x1 * y3 * z2 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1);
    }
    // 过三点求面方程通用式 Ax+By+Cx+D=0
    // 返回数组 [A,B,C,D]
    template <typename T = double, int M>
    btypes::Array<T> find_equation(
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2,
        const btypes::Array<T, M>& pt3)
    {
        assert(pt1.Size() == 3 && pt1.Size() == pt2.Size() && pt2.Size() == pt3.Size());
        return find_equation(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
    }
    // 过三点求面方程通用式 Ax+By+Cx+D=0
    // 以引用返回
    template <typename T = double, int M>
    void find_equation(
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2,
        const btypes::Array<T, M>& pt3,
        T& A, T& B, T& C, T& D)
    {
        assert(pt1.Size() == 3 && pt1.Size() == pt2.Size() && pt2.Size() == pt3.Size());
        btypes::Array<T> eq = find_equation(pt1, pt2, pt3);
        A = eq[0];
        B = eq[1];
        C = eq[2];
        D = eq[3];
    }
    // 求点到面的垂足, c1,c2,c3为面上任意三点, 不共线
    template <typename T = double, int M>
    btypes::Array<T> find_footpoint(
        const btypes::Array<T, M>& pt,
        const btypes::Array<T, M>& c1,
        const btypes::Array<T, M>& c2,
        const btypes::Array<T, M>& c3)
    {
        // 法向
        btypes::Array<T> n = cross(c2 - c1, c3 - c1);
        if (length(n) == 0)
            return btypes::Array<T>();
        T t = (n[0] * c1[0] + n[1] * c1[1] + n[2] * c1[2] - (n[0] * pt[0] + n[1] * pt[1] + n[2] * pt[2])) / (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        return btypes::Array<T>(pt[0] + n[0] * t, pt[1] + n[1] * t, pt[2] + n[2] * t);
    }
    /// \param 求这个点到平面的垂足
    /// \param A,B,C,D 平面方程 Ax+By+Cz+D = 0
    template <typename T = double, int M>
    btypes::Array<T> find_footpoint(
        const btypes::Array<T, M>& pt,
        const T& A, const T& B, const T& C, const T& D)
    {
        T t = (A * pt[0] + B * pt[1] + C * pt[2] + D) / (A * A + B * B + C * C);
        return btypes::Array<T>(
            pt[0] - A * t,
            pt[1] - B * t,
            pt[2] - C * t);
    }
    // 点到直线的垂足
    template <typename T = double, int M>
    btypes::Array<T> find_footpoint(
        const btypes::Array<T, M>& pt0,
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2)
    {
        assert(pt0.Size() == 3 && pt0.Size() == pt1.Size() && pt1.Size() == pt2.Size());
        double part1 = (btypes::Mat<T>((pt0 - pt1).Data(), 1, 3) * btypes::Mat<T>((pt2 - pt1).Data(), 3, 1))(0, 0);
        double part2 = (btypes::Mat<T>((pt2 - pt1).Data(), 1, 3) * btypes::Mat<T>((pt2 - pt1).Data(), 3, 1))(0, 0);
        if (part2 == 0)
            return btypes::Array<T>(); // 无解.
        else {
            double t = part1 / part2;
            return btypes::Array<T>(
                pt1[0] + (pt2[0] - pt1[0]) * t,
                pt1[1] + (pt2[1] - pt1[1]) * t,
                pt1[2] + (pt2[2] - pt1[2]) * t);
        }
    }
    // 过两点求直线方程, 参数式, 返回值为一个序列, [A,a, B,b, C,c]
    // x = A + a*t
    // y = B + b*t
    // z = C + c*t
    template <typename T = double, int M>
    btypes::Array<T> find_equation(
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2)
    {
        return btypes::Array<T>(pt1[0], pt2[0] - pt1[0], pt1[1], pt2[1] - pt1[1], pt1[2], pt2[2] - pt1[2]);
    }
    // 求直线交点
    template <typename T = double, int M>
    btypes::Array<T> find_crosspoint(
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2,
        const btypes::Array<T, M>& pt3,
        const btypes::Array<T, M>& pt4)
    {
        // x = A + a*t  = A' + a'*t
        // y = B + b*t  = B' + b'*t
        // z = C + c*t  = C' + c'*t
        // A + a*t + B + b*t + C + c*t = A' + a'*t + B' + b'*t + C' + c'*t
        // (a-a' + b-b' + c-c')*t = A'+B'+C'-A-B-C
        // t = (A'+B'+C'-A-B-C)/(a-a'+b-b'+c-c')
        btypes::Array<T> eq1 = find_equation(pt1, pt2);
        btypes::Array<T> eq2 = find_equation(pt3, pt4);
        double part1 = eq2[0] + eq2[2] + eq2[4] - eq1[0] - eq1[2] - eq2[4];
        double part2 = eq1[1] - eq2[1] + eq1[3] - eq2[3] + eq1[5] - eq2[5];
        if (part2 == 0)
            return btypes::Array<T>(); // 无解.
        else {
            double t = part1 / part2;
            return btypes::Array<T>(
                pt1[0] + (pt2[0] - pt1[0]) * t,
                pt1[1] + (pt2[1] - pt1[1]) * t,
                pt1[2] + (pt2[2] - pt1[2]) * t);
        }
    }
    // 求直线与面的交点
    // 求交点
    // 直线方程:
    //           (x-p1.x) / (p2.x - p1.x) = c
    //           (y-p1.y) / (p2.y - p1.y) = c
    //           (z-p1.z) / (p2.z - p1.z) = c
    //
    // 写成参数式
    // x = p1.x + (p2.x - p1.x) * t
    // y = p1.y + (p2.y - p1.y) * t
    // z = p1.z + (p2.z - p1.z) * t
    //
    // Ax + By + Cz + D = 0; 将x,y,z代入得到
    //
    // A(p1.x + (p2.x - p1.x)*t) + B(p1.y + (p2.y-1.y)*t) + C(p1.z + (p2.z-p1.z)*t) + D = 0
    // 求出 t 即可
    // A*p1.x + A*(p2.x-p1.x)*t + B....
    //
    // (A*(p2.x-p1.x) + B*(p2.y-p1.y) + C*(p2.z-p1.z))*t = -(A*p1.x + B*p1.y + C*p1.z + D)
    template <typename T = double, int M>
    btypes::Array<T> find_crosspoint(
        const btypes::Array<T, M>& pt1,
        const btypes::Array<T, M>& pt2,
        const T& A, const T& B, const T& C, const T& D)
    {
        double part1 = -(A * pt1[0] + B * pt1[1] + C * pt1[2] + D);
        double part2 = (A * (pt2[0] - pt1[0]) + B * (pt2[1] - pt1[1]) + C * (pt2[2] - pt1[2]));
        if (part2 == 0) {
            // 无解, 说明直线平行于平面
            return btypes::Array<T>();
        } else {
            double t = part1 / part2;
            return btypes::Array<T>(
                pt1[0] + (pt2[0] - pt1[0]) * t,
                pt1[1] + (pt2[1] - pt1[1]) * t,
                pt1[2] + (pt2[2] - pt1[2]) * t);
        }
    }
}
}

#endif
