///
/// (c) 2018 chivstyle, All rights reserved
///
#ifndef _nurbs_curve3d_h
#define _nurbs_curve3d_h

#include <stdio.h>

#ifdef _MSC_VER
#ifdef LIBNURBS_EXPORTS
#define LIBNURBS_API __declspec(dllexport)
#else
#define LIBNURBS_API __declspec(dllimport)
#endif
#else
#define LIBNURBS_API
#endif
/// \class NurbsCurve3D
///        这个是Nurbs的一个特例, B-spline, 首尾过控制点
///        1. 不支持动态增减控制点和节点向量, 但可动态更改
class LIBNURBS_API NurbsCurve3D {
public:
    NurbsCurve3D();
    NurbsCurve3D(const NurbsCurve3D&);
    NurbsCurve3D(NurbsCurve3D&&);
	NurbsCurve3D(int cv_count, int degree);
	virtual ~NurbsCurve3D();
    //
    bool IsValid() const;
    //
    void Create(int cv_count, int degree);
	// filename - Local encoding
	void DumpTo(const char *filename);
    void DumpTo(FILE* fp);
	//
	void SetKnot(int k, double knot);
	//
	void MakeClampedUniformKnotVector();
	// DO NOT USE IT
	void MakePeriodicUniformKnotVector();
	//
	void SetCv(int k, const double *p);
	void SetCv(int k, double x, double y, double z);
	//
	double *GetCv() { return mCv; }
	double *GetKnots() { return mKnots; }
    const double* GetCv() const { return mCv; }
    const double *GetKnots() const { return mKnots; }
	//
	int GetCvCount() const { return mCvCount; }
	int GetKnotCount() const { return mKnotCount; }
    int GetDegree() const { return mDegree; }
	// 生成曲线段
	void GenerateEv(double *ev, int point_count) const;
	// 从EV(型值点)生成BSpline.
	static NurbsCurve3D *CreateFromEv(const double *ev, int point_count, int degree);
    //
    static NurbsCurve3D FromEv(const double *ev, int point_count, int degree);
    //
    NurbsCurve3D& operator=(const NurbsCurve3D&);
    NurbsCurve3D& operator=(NurbsCurve3D&&);
	//
protected:
	// knot vector
	double *mKnots;
	// control points.
	double *mCv;
	//
	int mCvCount;
	//
	int mDegree;
	//
	int mKnotCount;
};

#endif
