///
/// (C) 2018 CHIV
///
/// @file btypes_t.h
/// @author chiv
/// @brief 提供一些方便使用的类, 如果需要使用专业的运动学或者线代运算, 这里的代码并不能够提供足够
///        的性能.
///
#pragma once

#include <assert.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <float.h>

#ifndef ARRAYSZ
#define ARRAYSZ(a) (sizeof(a) / sizeof(a[0]))
#endif
#ifndef B_PI
#define B_PI 3.1415926535897932384626
#endif

namespace btypes {

template <typename T>
static __inline T& _limit2(T& var, const T& lower, const T& upper)
{
    if (var < lower)
        var = lower;
    else if (var > upper)
        var = upper;
    return var;
}
template <typename T>
static __inline T D2R(T d)
{
    return static_cast<T>(d / 180. * B_PI);
}
template <typename T>
static __inline T R2D(T r)
{
    return static_cast<T>(r / B_PI * 180.);
}
template <typename T>
static inline T limit(T x, T lower, T upper)
{
    if (x < lower)
        x = lower;
    if (x > upper)
        x = upper;
    return x;
}
/// @param T typename
/// @param D Size of cache, unit : sizeof(T)
template <typename T, int D = 8>
class Array {
public:
    Array()
        : mArray(nullptr)
        , mCount(0)
    {
        mArray = Alloc(0);
    }
    //
    explicit Array(const T& _1)
        : mArray(nullptr)
        , mCount(1)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
    }
    //
    explicit Array(const T& _1, const T& _2)
        : mArray(nullptr)
        , mCount(2)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3)
        : mArray(nullptr)
        , mCount(3)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4)
        : mArray(nullptr)
        , mCount(4)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4, const T& _5)
        : mArray(nullptr)
        , mCount(5)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
        mArray[4] = _5;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4, const T& _5, const T& _6)
        : mArray(nullptr)
        , mCount(6)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
        mArray[4] = _5;
        mArray[5] = _6;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4, const T& _5, const T& _6, const T& _7)
        : mArray(nullptr)
        , mCount(7)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
        mArray[4] = _5;
        mArray[5] = _6;
        mArray[6] = _7;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4, const T& _5, const T& _6, const T& _7,
        const T& _8)
        : mArray(nullptr)
        , mCount(8)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
        mArray[4] = _5;
        mArray[5] = _6;
        mArray[6] = _7;
        mArray[7] = _8;
    }
    //
    explicit Array(const T& _1, const T& _2, const T& _3, const T& _4, const T& _5, const T& _6, const T& _7,
        const T& _8, const T& _9)
        : mArray(nullptr)
        , mCount(9)
    {
        mArray = Alloc(mCount);
        mArray[0] = _1;
        mArray[1] = _2;
        mArray[2] = _3;
        mArray[3] = _4;
        mArray[4] = _5;
        mArray[5] = _6;
        mArray[6] = _7;
        mArray[7] = _8;
        mArray[8] = _9;
    }
    //
    explicit Array(const T* arr, int count)
        : mCount(count)
        , mArray(nullptr)
    {
        mArray = Alloc(count);
        if (arr) {
            for (int k = 0; k < count; ++k) {
                mArray[k] = arr[k];
            }
        }
    }
    //
    template <int M>
    Array(const Array<T, M>& arr)
        : mCount(0)
        , mArray(0)
    {
        Free(mArray);
        mArray = Alloc(arr.Count());
        mCount = arr.Count();
        for (int k = 0; k < arr.Count(); ++k) {
            mArray[k] = arr[k];
        }
    }
    Array(const Array& arr)
        : mCount(0)
        , mArray(0)
    {
        Free(mArray);
        mArray = Alloc(arr.Count());
        mCount = arr.Count();
        for (int k = 0; k < arr.Count(); ++k) {
            mArray[k] = arr[k];
        }
    }
    //
    template <int M>
    Array(Array<T, M>&& arr)
        : mArray(0)
        , mCount(0)
    {
        Free(mArray);
        if (arr._Array() == arr._Cache()) {
            mArray = Alloc(arr.Count());
            for (int k = 0; k < arr.Count(); ++k) {
                mArray[k] = (T &&) arr[k];
            }
        } else {
            mArray = arr._Array();
        }
        mCount = arr.Count();
        mCapacity = arr.Capacity();
        arr._Count() = 0;
        arr._Array() = 0;
    }
    Array(Array&& arr)
        : mArray(0)
        , mCount(0)
    {
        Free(mArray);
        if (arr._Array() == arr._Cache()) {
            mArray = Alloc(arr.Count());
            for (int k = 0; k < arr.Count(); ++k) {
                mArray[k] = (T &&)arr[k];
            }
        } else {
            mArray = arr._Array();
        }
        mCount = arr.Count();
        mCapacity = arr.Capacity();
        arr._Count() = 0;
        arr._Array() = 0;
    }
    /*!
		* \brief resize the buffer
		* \warning This method will destroy your data.
		*/
    void Resize(int count)
    {
        Free(mArray);
        mCount = count;
        mArray = Alloc(mCount);
    }
    // 用于支持任意缓冲的array
    template <int M>
    void operator=(const Array<T, M>& arr)
    {
        Free(mArray);
        mArray = Alloc(arr.Count());
        mCount = arr.Count();
        for (int k = 0; k < arr.Count(); ++k) {
            mArray[k] = arr.Data()[k];
        }
    }
    void operator=(const Array& arr)
    {
        Free(mArray);
        mArray = Alloc(arr.Count());
        mCount = arr.Count();
        for (int k = 0; k < arr.Count(); ++k) {
            mArray[k] = arr.Data()[k];
        }
    }
    //
    template <int M>
    void operator=(Array<T, M>&& arr)
    {
        Free(mArray);
        if (arr._Array() == arr._Cache()) {
            mArray = Alloc(arr.Count());
            for (int k = 0; k < arr.Count(); ++k) {
                mArray[k] = (T &&) arr.Data()[k];
            }
        } else {
            mArray = arr._Array();
        }
        mCount = arr.Count();
        mCapacity = arr.Capacity();
        arr._Count() = 0;
        arr._Array() = 0;
    }
    void operator=(Array&& arr)
    {
        Free(mArray);
        if (arr._Array() == arr._Cache()) {
            mArray = Alloc(arr.Count());
            for (int k = 0; k < arr.Count(); ++k) {
                mArray[k] = (T &&)arr.Data()[k];
            }
        } else {
            mArray = arr._Array();
        }
        mCount = arr.Count();
        mCapacity = arr.Capacity();
        arr._Count() = 0;
        arr._Array() = 0;
    }
    //
    template <int M>
    Array operator-(const Array<T, M>& arr) const
    {
        assert(Count() == arr.Count());
        Array<T> t;
        t.Resize(arr.Count());
        for (int k = 0; k < t.Count(); ++k) {
            t[k] = mArray[k] - arr[k];
        }
        return (Array &&) t;
    }
    //
    template <int M>
    Array operator+(const Array<T, M>& arr) const
    {
        assert(Count() == arr.Count());
        Array<T> t;
        t.Resize(arr.Count());
        for (int k = 0; k < t.Count(); ++k) {
            t[k] = mArray[k] + arr[k];
        }
        return (Array &&) t;
    }
    Array operator/(const T& d) const
    {
        Array<T> t;
        t.Resize(Count());
        for (int k = 0; k < t.Count(); ++k) {
            t[k] = mArray[k] / d;
        }
        return (Array &&) t;
    }
    //
    Array Slice(int from, int len = -1)
    {
        if (len == -1) {
            len = mCount - from;
        }
        btypes::Array slice; slice.Resize(len);
        for (int k = 0; k < len; ++k) {
            slice[k] = mArray[from + k];
        }
        return (Array&&)slice;
    }
    //
    Array& operator<<(const T& d)
    {
        return Push(d);
    }
    Array& operator,(const T& d)
    {
        return Push(d);
    }
    template <int M>
    Array& operator<<(const Array<T, M>& d)
    {
        return Push(d);
    }
    template <int M>
    Array& operator,(const Array<T, M>& d)
    {
        return Push(d);
    }
    //
    template <int M>
    Array& Push(const Array<T, M>& a)
    {
        if (mCapacity >= mCount + a.Count()) {
            for (int k = 0; k < a.Count(); ++k) {
                mArray[mCount++] = a[k];
            }
        } else {
            T* _old_buffer = mArray;
            mArray = Alloc(mCount + a.Count());
            for (int k = 0; k < mCount; ++k) {
                mArray[k] = _old_buffer[k];
            }
            for (int k = 0; k < a.Count(); ++k) {
                mArray[mCount + k] = a[k];
            }
            mCount += a.Count();
            Free(_old_buffer);
        }
        return *this;
    }
    //
    Array& Push(const T& d)
    {
        if (mCapacity > mCount) { // >= mCount+1
            mArray[mCount++] = d;
        } else {
            // we need more memory to store data.
            T* _old_buffer = mArray;
            mArray = Alloc(mCount+1);
            for (int k = 0; k < mCount; ++k) {
                mArray[k] = _old_buffer[k];
            }
            mArray[mCount++] = d;
            Free(_old_buffer);
        }
        return *this;
    }
    //
    virtual ~Array()
    {
        Free(mArray);
    }
    // clear allocated resource
    // release_memory - If set this parameter true, we'll clear
    // the data then release the memory allocated.
    void Clear(bool release_memory = false)
    {
        if (release_memory) {
            Free(mArray); mArray = nullptr;
            mCapacity = 0;
        }
        mCount = 0;
    }
    //
    int Capacity() const
    {
        return mCapacity;
    }
    //
    int Count() const
    {
        return mCount;
    }
    //
    int Size() const
    {
        return mCount;
    }
    Array operator*(const T& t)
    {
        Array<T> ba;
        ba.Resize(Size());
        for (int k = 0; k < Size(); ++k) {
            ba[k] = mArray[k] * t;
        }
        // move
        return (Array &&) ba;
    }
    //
    T& operator[](int index)
    {
        assert(index >= 0 && index < mCount && "index is not valid");
        return mArray[index];
    }
    //
    const T& operator[](int index) const
    {
        assert(index >= 0 && index < mCount && "index is not valid");
        return mArray[index];
    }
    //
    const T* Data() const
    {
        return mArray;
    }
    //
    T* Data()
    {
        return mArray;
    }
    // 不要从外部调用它们.
    int& _Count() { return mCount; }
    T*& _Array() { return mArray; }
    T* _Cache() { return mBuffer; }
    //
protected:
    //
    T* Alloc(int count)
    {
        if (count > ARRAYSZ(mBuffer)) { // Cache is smaller than count
                                        // we should allocate memory on heap.
            mCapacity = 16 + count;     // 16 more than requested count.
            // Allocate more.
            return new T[mCapacity];
        } else {
            mCapacity = ARRAYSZ(mBuffer);
            return &mBuffer[0];
        }
    }
    //
    void Free(T* t)
    {
        if (t != &mBuffer[0] && t) {
            delete[] t;
        }
    }

private:
    //
    int mCount;
    int mCapacity;
    // buffer.
    T* mArray;
    // chache.
    T mBuffer[D];
};
// 行主序
template <typename T, int D = 16>
class Mat {
public:
    Mat()
        : mShiftCount(0)
    {
        Resize(0, 0);
    }
    //
    Mat(int row_count, int col_count)
        : mShiftCount(0)
    {
        Resize(row_count, col_count);
    }
    //
    Mat(const T* data, int row_count, int col_count)
        : mArray(data, col_count*row_count)
        , mShiftCount(0)
    {
        mCount[0] = row_count;
        mCount[1] = col_count;
    }
    //
    void Resize(int row_count = 3, int col_count = 3)
    {
        mArray.Resize(col_count * row_count);
        mCount[0] = row_count;
        mCount[1] = col_count;
    }
    //
    void SetIdentity()
    {
        for (int k = 0; k < mCount[0]; ++k) {
            for (int u = 0; u < mCount[1]; ++u) {
                mArray[mCount[1] * k + u] = u == k ? 1 : 0;
            }
        }
    }
    //
    void SetZero()
    {
        for (int k = 0; k < mCount[0]; ++k) {
            for (int u = 0; u < mCount[1]; ++u) {
                mArray[mCount[1] * k + u] = 0;
            }
        }
    }
    //
    void UpdateRow(int row_, const Array<T>& data)
    {
        assert(mCount[1] == data.Count() && row_ < mCount[0]);
        //
        for (int k = 0; k < data.Count(); ++k) {
            mArray[mCount[1] * row_ + k] = data[k];
        }
    }
    //
    void UpdateCol(int col_, const Array<T>& data)
    {
        assert(mCount[0] == data.Count() && col_ < mCount[1]);
        //
        for (int k = 0; k < data.Count(); ++k) {
            mArray[mCount[1] * k + col_] = data[k];
        }
    }
    //
    btypes::Array<double> GetCol(int n) const
    {
        assert(n >= 0 && n < mCount[1]);
        btypes::Array<double> column;
        column.Resize(mCount[0]);
        for (int k = 0; k < mCount[0]; ++k) {
            column[k] = mArray[mCount[1] * k + n];
        }
        return column;
    }
    btypes::Array<double> GetRow(int n) const
    {
        assert(n >= 0 && n < mCount[0]);
        btypes::Array<double> row;
        row.Resize(mCount[1]);
        for (int k = 0; k < mCount[1]; ++k) {
            row[k] = mArray[mCount[1] * n + k];
        }
        return row;
    }
    //
    T& operator()(int row_, int col_)
    {
        assert(row_ < Count()[0] && col_ < Count()[1]);
        //
        return mArray[mCount[1] * row_ + col_];
    }
    //
    const T& operator()(int row_, int col_) const
    {
        assert(row_ < Count()[0] && col_ < Count()[1]);
        //
        return mArray[mCount[1] * row_ + col_];
    }
    //
    Mat& operator<<(const T& d)
    {
        assert(mShiftCount < mArray.Size());
        //
        mArray[mShiftCount++] = d;
        return *this;
    }
    Mat& operator,(const T& d)
    {
        assert(mShiftCount < mArray.Size());
        //
        mArray[mShiftCount++] = d;
        return *this;
    }
    //
    static Mat Transpose(const Mat& mat_)
    {
        Mat mat(mat_.Count()[1], mat_.Count()[0]);
        for (int k = 0; k < mat_.Count()[0]; ++k) {
            for (int u = 0; u < mat_.Count()[1]; ++u) {
                mat(u, k) = mat_(k, u);
            }
        }
        return (Mat &&) mat;
    }
    //
    void operator*=(const T& d)
    {
        for (int k = 0; k < mCount[0]; ++k) {
            for (int u = 0; u < mCount[1]; ++u) {
                (*this)(k, u) *= d;
            }
        }
    }
    Mat operator*(const T& d) const
    {
        Mat mat_ = *this;
        for (int k = 0; k < mCount[0]; ++k) {
            for (int u = 0; u < mCount[1]; ++u) {
                mat_(k, u) *= d;
            }
        }
        return (Mat &&) mat_;
    }
    //
    template <int M>
    Mat operator-(const Mat<T, M>& mat_) const
    {
        assert(Count()[0] == mat_.Count()[0] && Count()[1] == mat_.Count()[1]);
        //
        return Mat((mArray - mat_._Array()).Data(), Count()[0], Count()[1]);
    }
    template <int M>
    Mat operator+(const Mat<T, M>& mat_) const
    {
        assert(Count()[0] == mat_.Count()[0] && Count()[1] == mat_.Count()[1]);
        //
        return Mat((mArray + mat_._Array()).Data(), Count()[0], Count()[1]);
    }
    //
    template <int M>
    Mat operator*(const Mat<T, M>& mat_) const
    {
        assert(mCount[1] == mat_.mCount[0]);
        Mat result(mCount[0], mat_.Count(1));
        result.SetZero(); // 设置0
        // let's go.
        // 遍历结果行, 一个一个求结果
        for (int u = 0; u < mCount[0]; ++u) {
            for (int v = 0; v < mat_.Count(1); ++v) {
                for (int k = 0; k < mCount[1]; ++k) {
                    // u行 乘 v列
                    result(u, v) += (*this)(u, k) * mat_(k, v);
                }
            }
        }
        //
        return (Mat &&)result;
    }
    enum MatrixError {
        E_OK,
        E_BAD_MATRIX,
        E_NO_RESOURCE,
        E_NO_INV
    };
    static int matrix_inverse(const double *dmatrix, int row, int col, double *out)
    {
        /** temp vars */
        int    *is, *js, i, j, k, l, u, v, n;
        double  d, p;
        double *matrix;
        /** check the matrix, must be square */
        if (row != col) return -E_BAD_MATRIX;
        n = row;
        is = (int *)malloc(sizeof(int) * row);
        if (!is) return -E_NO_RESOURCE;
        js = (int *)malloc(sizeof(int) * row);
        if (!js) { free(is); return -E_NO_RESOURCE; }
        matrix = (double *)malloc(sizeof(double) * n * n);
        if (!matrix) { free(is);free(js);return -E_NO_RESOURCE; }
        memcpy(matrix, dmatrix, sizeof(double) * n * n);
        /** okay, resources ready */
        for (k = 0; k < n; k++) {
            d = 0.;
            for (i = k; i < n; i++) {
                for (j = k; j < n; j++) {
                    l = i * n + j;
                    p = fabs(matrix[l]);
                    if (p > d) {
                        d = p;
                        is[k] = i;
                        js[k] = j;
                    }
                }
            }
            if (d == 0.) {
                free(is);
                free(js);
#if defined(_DEBUG)
                memcpy(out, matrix, sizeof(double) * n * n);
#endif
                free(matrix);
                return -E_NO_INV;
            }
            if (is[k] != k) {
                for (j = 0; j < n; j++) {
                    u = k * n + j;
                    v = is[k] * n + j;
                    p = matrix[u];
                    matrix[u] = matrix[v];
                    matrix[v] = p;
                }
            }
            if (js[k] != k) {
                for (i = 0; i < n; i++) {
                    u = i * n + k;
                    v = i * n + js[k];
                    p = matrix[u];
                    matrix[u] = matrix[v];
                    matrix[v] = p;
                }
            }
            l = k * n + k;
            matrix[l] = 1. / matrix[l];
            for (j = 0; j < n; j++) {
                if (j != k) {
                    u = k * n + j;
                    matrix[u] = matrix[u] * matrix[l];
                }
            }
            for (i = 0; i < n; i++) {
                if (i != k) {
                    for (j = 0; j < n; j++) {
                        if (j != k) {
                            u = i * n + j;
                            matrix[u] = matrix[u] - matrix[i * n + k] * matrix[k * n + j];
                        }
                    }
                    u = i * n + k;
                    matrix[u] = -matrix[u] * matrix[l];
                }
            }
        }
        for (k = n - 1; k >= 0; k--) {
            if (js[k] != k) {
            for (j = 0; j < n; j++) {
                    u = k * n + j;
                    v = js[k] * n + j;
                    p = matrix[u];
                    matrix[u] = matrix[v];
                    matrix[v] = p;
                }
            }
            if (is[k] != k) {
                for (i = 0; i < n; i++) {
                    u = i * n + k;
                    v = i * n + is[k];
                    p = matrix[u];
                    matrix[u] = matrix[v];
                    matrix[v] = p;
                }
            }
        }
        memcpy(out, matrix, sizeof(double) * n * n);
        free(is);
        free(js);
        free(matrix);

        return -E_OK;
    }
    /// from internet, LUdecompose
    /// 虽然有些代码与本程序中的某部分功能重复, 但我们尊重原始程序, 只更改风格, 名称
    /// 其余的保留原样, 这样做的目的只是为了风格统一, 并无他意. 另外加了奇异判断.
    Mat<double, D> Inverse() const
    {
        Mat<double, D> inv(Count()[0], Count()[1]);
        if (matrix_inverse(Data(), Count()[0], Count()[1], inv.Data()) != E_OK) {
            return Mat<double, D>(); // Empty, stands for failure.
        }
        return inv;
    }
    //
    Mat Transpose()
    {
        *this = (Mat&&)Transpose(*this);
        //
        return *this;
    }
    btypes::Array<int, 2> Size() const
    {
        return btypes::Array<int, 2>(mCount[0], mCount[1]);
    }
    // Convenient methods.
    __inline int Count(int i) const { assert(i >= 0 && i < sizeof(mCount)); return mCount[i]; }
    __inline int NumOfCols() const { return mCount[1]; }
    __inline int NumOfRows() const { return mCount[0]; }
    //
    __inline const int* Count() const { return mCount; }
    //
    __inline const T* Data() const { return mArray.Data(); }
    __inline T* Data() { return mArray.Data(); }
    //
    __inline const Array<T, D>& _Array() const { return mArray; }
    __inline Array<T, D>& _Array() { return mArray; }
    //
private:
    Array<T, D> mArray;
    //
    int mCount[2]; // 0-row, 1-col
    int mShiftCount;
};
//
struct POSE {
    double x, y, z;
    double Rx, Ry, Rz;
};
//
static __inline Mat<double> Make(const struct POSE& pose)
{
    double c[] = { cos(D2R(pose.Rx)), cos(D2R(pose.Ry)), cos(D2R(pose.Rz)) };
    double s[] = { sin(D2R(pose.Rx)), sin(D2R(pose.Ry)), sin(D2R(pose.Rz)) };

    Mat<double> transform(4, 4);
    transform.UpdateRow(0, Array<double>(c[1] * c[2], s[0] * s[1] * c[2] - c[0] * s[2], c[0] * s[1] * c[2] + s[0] * s[2], pose.x));
    transform.UpdateRow(1, Array<double>(c[1] * s[2], s[0] * s[1] * s[2] + c[0] * c[2], c[0] * s[1] * s[2] - s[0] * c[2], pose.y));
    transform.UpdateRow(2, Array<double>(-s[1], s[0] * c[1], c[0] * c[1], pose.z));
    transform.UpdateRow(3, Array<double>(0, 0, 0, 1));
    //
    return (Mat<double> &&) transform;
}
//
static __inline Mat<double> Make(double x, double y, double z, double rx, double ry, double rz)
{
    return Make(POSE{ x, y, z, rx, ry, rz });
}
//
static __inline POSE Make(const Mat<double>& mat_)
{
    const double &R31 = mat_(2, 0), &R32 = mat_(2, 1), &R33 = mat_(2, 2),
	             &R21 = mat_(1, 0),
                 &R11 = mat_(0, 0), &R12 = mat_(0, 1), &R13 = mat_(0, 2);
    double theta = -asin(btypes::limit(R31, -1., 1.));
    double cos_theta = cos(theta);
    double psi = 0, phi = 0.;
    if (fabs(cos_theta) < FLT_MIN) {
        phi = 0.; // any value will be ok.
        if (theta > 0) { // Pi/2
            psi = phi + atan2(R12, R13);
        } else { // -Pi/2
            psi = -phi + atan2(-R12, -R13);
        }
    } else {
        psi = atan2(R32 / cos_theta, R33 / cos_theta);
        phi = atan2(R21 / cos_theta, R11 / cos_theta);
    }
    //
    return POSE{ mat_(0, 3), mat_(1, 3), mat_(2, 3), 
        R2D(psi), R2D(theta), R2D(phi) };
}
// Quat
class Quat {
public:
    Quat(double x, double y, double z, double w)
    {
        mQuat[0] = x;
        mQuat[1] = y;
        mQuat[2] = z;
        mQuat[3] = w;
    }
    //
    Quat(double rx /*psi*/, double ry /*theta*/, double rz /*phi*/)
    {
        double c[] = {
            cos(D2R(rx) / 2), cos(D2R(ry) / 2), cos(D2R(rz) / 2)
        };
        double s[] = {
            sin(D2R(rx) / 2), sin(D2R(ry) / 2), sin(D2R(rz) / 2)
        };
        mQuat[0] = c[2] * c[1] * c[0] + s[2] * s[1] * s[0];
        mQuat[1] = s[2] * c[1] * c[0] - c[2] * s[1] * s[0];
        mQuat[2] = c[2] * s[1] * c[0] + s[2] * c[1] * s[0];
        mQuat[3] = c[2] * s[1] * s[0] - s[2] * s[1] * c[0];
    }
    //
    Quat(const Array<double>& v, double angle)
    {
        double rad = D2R(angle);
        double c = cos(rad / 2);
        double s = sin(rad / 2);
        mQuat[3] = c;
        mQuat[0] = s * v[0];
        mQuat[1] = s * v[1];
        mQuat[2] = s * v[2];
    }
    //
    Array<double> ToEuler()
    {
        const double& x_ = mQuat[0];
        const double& y_ = mQuat[1];
        const double& z_ = mQuat[2];
        const double& w_ = mQuat[3];
        //
        double psi = atan2(2 * (w_ * x_ + y_ * z_), 1 - 2 * (x_ * x_ + y_ * y_));
        double theta = asin(limit<double>(2 * (w_ * y_ - z_ * x_), -1, 1));
        double phi = atan2(2 * (w_ * z_ + x_ * y_), 1 - 2 * (y_ * y_ + z_ * z_));
        //
        return Array<double>(R2D(psi), R2D(theta), R2D(phi));
    }
    //
    double x() const { return mQuat[0]; }
    double y() const { return mQuat[1]; }
    double z() const { return mQuat[2]; }
    double w() const { return mQuat[3]; }
    //
    double mQuat[4];
};
} // namespace BTypes
