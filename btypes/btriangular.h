///
/// (C) 2018 chiv, All rights reserved.
///
/// \brief Need C++11 support.
///
#pragma once

#include <cassert>
#include <list>
#include <vector>

namespace btypes {

    namespace btriangular {
        ///
        enum ScanDirection {
            CLOCKWISE,
            ANTICLOCKWISE,
            BOOMWISE // WTF.
        };
        template <typename T>
        class Vector3d {
        public:
            Vector3d() {}
            Vector3d(T x_, T y_, T z_)
                : x(x_)
                , y(y_)
                , z(z_)
            {
            }
            T x, y, z;
        };
        template <typename T>
        class Point2d {
        public:
            Point2d() {}
            Point2d(T x_, T y_)
                : x(x_)
                , y(y_)
            {
            }
            T x, y;
        };
        template <typename T>
        class Point3d {
        public:
            Point3d() {}
            Point3d(T x_, T y_, T z_)
                : x(x_)
                , y(y_)
                , z(z_)
            {
            }
            T x, y, z;
        };
        template <typename DataType, typename IndexType = int>
        class Point3di : public Point3d<DataType> {
        public:
            Point3di() {}
            Point3di(DataType x_, DataType y_, DataType z_, IndexType idx_)
                : Point3d<DataType>(x_, y_, z_)
                , idx(idx_)
            {
            }
            IndexType idx;
        };
        ///
        template <typename T>
        Vector3d<T> CrossProduct(const Vector3d<T>& v1, const Vector3d<T>& v2)
        {
            Vector3d<T> pd = {
                v1.y * v2.z - v2.y * v1.z, // be 0 in plane2d
                v1.z * v2.x - v1.x * v2.z, // be 0 in plane2d
                v1.x * v2.y - v2.x * v1.y
            };
            return pd;
        }
        ///
        template <typename T>
        bool IsPointInTriangles(const Point3d<T>& P, const Point3d<T>& A, const Point3d<T>& B, const Point3d<T>& C)
        {
            Vector3d<T> BC = { C.x - B.x, C.y - B.y, 0 }; // BC
            Vector3d<T> BP = { P.x - B.x, P.y - B.y, 0 }; // BP
            Vector3d<T> CA = { A.x - C.x, A.y - C.y, 0 }; // CA
            Vector3d<T> CP = { P.x - C.x, P.y - C.y, 0 }; // CP
            Vector3d<T> AB = { B.x - A.x, B.y - A.y, 0 }; // AB
            Vector3d<T> AP = { P.x - A.x, P.y - A.y, 0 }; // AP
            auto pp1 = CrossProduct(BC, BP);
            auto pp2 = CrossProduct(CA, CP);
            auto pp3 = CrossProduct(AB, AP);
            return (pp1.z <= 0 && pp2.z <= 0 && pp3.z <= 0) || (pp1.z >= 0 && pp2.z >= 0 && pp3.z >= 0);
        }
        /// \brief 将输入点组合成三角形簇
        /// \param points
        /// \param dims dimensions 可以是2或者3, 这个算法只针对平面, 因此如果输入了三维点
        ///        那么第三维将被忽略, dims的作用可以减少因为维度造成的缓冲重排.
        template <typename T>
        static __inline void _OUTPUT3D(const Point3d<T>& pt3d, std::vector<T>& output, int dims)
        {
            output.push_back(pt3d.x);
            output.push_back(pt3d.y);
            if (dims >= 3)
                output.push_back(pt3d.z);
        }

        template <typename T, int D = 3>
        ScanDirection GuessScanDirection(const T* points, int count, bool closed = false)
        {
            static_assert(D >= 2, __FUNCTION__ ":D should >= 2");
            if (count < 3)
                return BOOMWISE;
            if (closed)
                count = count - 1;
            int i0 = 0;
            for (int k = 1; k < count; ++k) {
                if (points[D * i0 + 0] > points[D * k + 0])
                    i0 = k;
            }
            int p_prev, p_next;
            if (i0 == 0) {
                p_prev = count - 1;
                p_next = 1;
            } else if (i0 == count - 1) {
                p_prev = count - 2;
                p_next = 0;
            } else {
                p_prev = i0 - 1;
                p_next = i0 + 1;
            }
            if (points[D * p_prev + 1] > points[D * p_next + 1]) {
                return ANTICLOCKWISE;
            } else if (points[D * p_prev + 1] < points[D * p_next + 1]) {
                return CLOCKWISE;
            } else {
                return BOOMWISE;
            }
        }
        /// \param points - 点缓冲, 一个数组, 如果是装箱类型, 要提供+,-等操作符.
        /// \param count  - 缓冲中点的数量, **是点的数量**
        /// \param closed - 是否是闭合的
        /// \param dims   - 点的维度, 必须>=2, 输出只会抽取有效点数, 如果dims是2, 则输出中每个点都是2维的
        ///                 如果dims>=3, 则每个点都是三维的, 其余的维度被忽略.
        /// \param direction 输入点的排布方向, 可以用GuessScanDirection获取.
        /// \return       - 返回一个数组, 包含重排后的三角形索引值.
        template <typename DataType, int D = 3, typename IndexType = int>
        std::vector<IndexType> Triangulate2(const DataType* points, int count, bool closed = false, int direction = ANTICLOCKWISE)
        {
            static_assert(D >= 2, __FUNCTION__ ":D should >= 2");
            //
            std::vector<IndexType> output;
            // push all points to list.
            std::vector<Point3di<DataType, IndexType>> vertices;
            if (count < 3)
                return output;
            for (int k = 0; k < count; ++k) {
                vertices.push_back(Point3di<DataType, IndexType>{
                    points[D * k + 0], points[D * k + 1], 0, k });
            }
            if (!closed) { // 如果曲线未闭合, 作三角形分割是没有意义的, 所以如果未闭合的话,
                // 会在计算时自动闭合.
                vertices.push_back(*vertices.begin());
            }
            // 把第二点加入队尾, 是为了形成 N,0,1的循环结构, 这就造成实际顶点数比真实的要多2个
            // 因此判定结束的条件就是顶点数==5, 其中前三个即最后一个三角形的三个顶点.
            vertices.push_back(*(vertices.begin() + 1));
            // let terror begin.
            do {
                // 每次都从新的输入的头部开始扫描
                std::vector<Point3di<DataType, IndexType>>::iterator it_s = vertices.begin();
                // find one.
                for (; it_s < vertices.end() - 2; ++it_s) {
                    const Point3di<DataType>&A = *(it_s + 0), &B = *(it_s + 1), &C = *(it_s + 2);
                    // if ABC valid ?
                    DataType nz = CrossProduct(Vector3d<DataType>{ B.x - A.x, B.y - A.y, B.z - A.z } /*AB*/,
                                               Vector3d<DataType>{ C.x - A.x, C.y - A.y, C.z - A.z } /*AC*/).z;
                    if (direction == ANTICLOCKWISE && nz <= 0 || direction == CLOCKWISE && nz >= 0)
                        continue;
                    // 别的点是否落在ABC构成的三角形内.
                    bool trapped = false;
                    // 前半拉
                    for (auto it = vertices.begin(); it < it_s; ++it) {
                        if (IsPointInTriangles(*it, A, B, C)) {
                            trapped = true;
                            break;
                        }
                    }
                    // 后半剌
                    if (!trapped) {
                        for (auto it = it_s + 3; it < vertices.end() - 2; ++it) {
                            if (IsPointInTriangles(*it, A, B, C)) {
                                trapped = true;
                                break;
                            }
                        }
                    }
                    if (trapped) continue; else break;
                }
                if (it_s == vertices.end() - 2) {
                    // 当线段发生交错时, 确实存在找不到凸角的情况, 此时直接返回一个错误的结果.
                    return std::move(output);
                }
                // 输出
                output.push_back((it_s + 0)->idx);
                output.push_back((it_s + 1)->idx);
                output.push_back((it_s + 2)->idx);
                // 移除凸点
                vertices.erase(it_s + 1);
            } while (vertices.size() > 5);
            output.push_back((vertices.begin() + 0)->idx);
            output.push_back((vertices.begin() + 1)->idx);
            output.push_back((vertices.begin() + 2)->idx);
            //
            return std::move(output);
        }
        //
        template <typename T, int D = 3>
        std::vector<T> Triangulate1(const T* points, int count, bool closed, int direction)
        {
            static_assert(D >= 2, __FUNCTION__ ":D should >= 2");
            std::vector<T> output;
            std::vector<int> indices = Triangulate2<T, D, int>(points, count, closed, direction);
            for each(auto idx in indices)
            {
                for (int k = 0; k < D; ++k) {
                    output.push_back(points[D * idx + k]);
                }
            }
            return std::move(output);
        }
        //
        template <typename T>
        __inline int _compare(T x1, T x2, T delta)
        {
            T dx = x1 - x2;
            if (dx > delta)
                return 1;
            else if (dx < delta)
                return -1;
            else
                return 0;
        }
        /// \brief 逆时针重排
        /// \param ev - 输入的点缓冲, 这里的点有序的闭合的曲线
        /// \param count
        /// \param out
        /// \param dims
        /// \param dx - 浮点比较精度
        template <typename T = double, int D = 2>
        int Arrange(const T* ev, int count, T* out, T dx = (T)(1e-10))
        {
            static_assert(D >= 2, __FUNCTION__ ":D should >= 2");
            //
            union _Pointx {
                T xd[D];
                struct {
                    T x, y;
                };
            }*p = (_Pointx*)ev, *o = (_Pointx*)out;
            // 找到最左边的点最靠上的点
            int i0 = 0;
            for (int k = 1; k < count; ++k) {
                // 比较的粗糙一点以获取更多的重复点
                int _r = _compare(p[i0].x, p[k].x, dx);
                if (_r > 0) {
                    i0 = k;
                } else if (_r == 0) {
                    if (_compare(p[i0].y, p[k].y, dx) < 0)
                        i0 = k;
                }
            }
            // 取数的方向
            int p_prev, p_next;
            if (i0 == 0) {
                p_prev = count - 1;
                p_next = 1;
            } else if (i0 == count - 1) {
                p_prev = count - 2;
                p_next = 0;
            } else {
                p_prev = i0 - 1;
                p_next = i0 + 1;
            }
            int op = 0;
            if (_compare(p[p_prev].y, p[p_next].y, dx) < 0) {
                // 朝前取
                for (int k = i0; k >= 0; --k) {
                    o[op++] = p[k];
                }
                for (int k = count - 1; k > i0; --k) {
                    o[op++] = p[k];
                }
            } else {
                // 朝后取
                for (int k = i0; k < count; ++k) {
                    o[op++] = p[k];
                }
                for (int k = 0; k < i0; ++k) {
                    o[op++] = p[k];
                }
            }
            return op;
        }
    } // btriangular

} // BTypes
