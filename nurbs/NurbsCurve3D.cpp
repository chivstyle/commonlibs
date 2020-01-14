///
/// (c) 2018 chiv, All rights reserved
///
#include "NurbsCurve3D.h"
#include "extra/NURBS.H"
#include "extra/interpolation.h"
#include "extra/matrix_error.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

NurbsCurve3D::NurbsCurve3D()
    : mCvCount(0)
    , mDegree(0)
    , mKnotCount(0)
    , mCv(nullptr)
    , mKnots(nullptr)
{
}

NurbsCurve3D::NurbsCurve3D(int cv_count, int degree)
    : NurbsCurve3D()
{
    Create(cv_count, degree);
}

void NurbsCurve3D::Create(int cv_count, int degree)
{
    if (mCvCount != cv_count || mDegree != degree) {
        delete[] mCv;
        delete[] mKnots;
        //
        mCvCount = cv_count;
        mDegree = degree;
        mKnotCount = cv_count + degree + 1;
        //
        mCv = new double[3 * mCvCount];
        mKnots = new double[mKnotCount];
    }
}

NurbsCurve3D::~NurbsCurve3D()
{
	delete[] mCv;
	delete[] mKnots;
}
NurbsCurve3D::NurbsCurve3D(const NurbsCurve3D& nc)
{
    mCvCount = nc.GetCvCount();
    mKnotCount = nc.GetKnotCount();
    mDegree = nc.GetDegree();
    if (nc.mCv && nc.mKnots) {
        mCv = new double[3 * nc.GetCvCount()];
        mKnots = new double[nc.GetKnotCount()];
        memcpy(mCv, nc.mCv, sizeof(double) * 3 * mCvCount);
        memcpy(mKnots, nc.mKnots, sizeof(double) * mKnotCount);
    }
}
NurbsCurve3D::NurbsCurve3D(NurbsCurve3D&& nc)
{
    mCvCount = nc.GetCvCount();
    mKnotCount = nc.GetKnotCount();
    mDegree = nc.GetDegree();
    mCv = nc.mCv;
    mKnots = nc.mKnots;
    //
    nc.mCv = nullptr;
    nc.mKnots = nullptr;
    nc.mCvCount = 0;
    nc.mKnotCount = 0;
    nc.mDegree = 0;
}
NurbsCurve3D& NurbsCurve3D::operator = (const NurbsCurve3D& nc)
{
    if (mCvCount != nc.GetCvCount() || mDegree != nc.GetDegree() ||
        mKnotCount != nc.GetKnotCount()) {
        delete[] mCv;
        delete[] mKnots;
        // allocate new.
        mCv = new double[3 * nc.GetCvCount()];
        mKnots = new double[nc.GetKnotCount()];
        mCvCount = nc.GetCvCount();
        mKnotCount = nc.GetKnotCount();
        mDegree = nc.GetDegree();
    }
    if (nc.mCv && nc.mKnots) {
        memcpy(mCv, nc.mCv, sizeof(double) * 3 * mCvCount);
        memcpy(mKnots, nc.mKnots, sizeof(double) * mKnotCount);
    }
    //
    return *this;
}
NurbsCurve3D& NurbsCurve3D::operator = (NurbsCurve3D&& nc)
{
    delete[] mCv;
    delete[] mKnots;
    //
    mCvCount = nc.GetCvCount();
    mKnotCount = nc.GetKnotCount();
    mDegree = nc.GetDegree();
    mCv = nc.mCv;
    mKnots = nc.mKnots;
    //
    nc.mCv = nullptr;
    nc.mKnots = nullptr;
    nc.mCvCount = 0;
    nc.mKnotCount = 0;
    nc.mDegree = 0;
    //
    return *this;
}

bool NurbsCurve3D::IsValid() const
{
    return mCvCount + mDegree + 1 == mKnotCount && mKnots && mCv;
}

void NurbsCurve3D::DumpTo(const char *filename)
{
    FILE* fp;
    errno_t es = fopen_s(&fp, filename, "wb");
    if (es == 0) {
        DumpTo(fp);
        fclose(fp);
    }
}

void NurbsCurve3D::DumpTo(FILE* fp)
{
    fprintf(fp, "Degree:%d\n", mDegree);
    fprintf(fp, "Cv:\n");
    for (int k = 0; k < mCvCount; ++k) {
        fprintf(fp, "%2d: [%16f, %16f, %16f]\n", k, mCv[3 * k + 0], mCv[3 * k + 1], mCv[3 * k + 2]);
    }
    fprintf(fp, "Knots:\n");
    for (int k = 0; k < mKnotCount; ++k) {
        fprintf(fp, "%2d: [%16f]\n", k, mKnots[k]);
    }
}

void NurbsCurve3D::MakeClampedUniformKnotVector()
{
	::MakeClampedUniformKnotsVector(mKnots, mKnotCount, mCvCount, mDegree);
}

void NurbsCurve3D::MakePeriodicUniformKnotVector()
{
	::MakePeriodicUniformKnotsVector(mKnots, mKnotCount, mCvCount, mDegree);
}

void NurbsCurve3D::SetKnot(int k, double knot)
{
	assert(k >= 0 && k < mKnotCount);
	//
	mKnots[k] = knot;
}

void NurbsCurve3D::SetCv(int k, const double *p)
{
	memcpy(&mCv[3 * k], p, sizeof(double) * 3);
}

void NurbsCurve3D::SetCv(int k, double x, double y, double z)
{
	mCv[3 * k + 0] = x;
	mCv[3 * k + 1] = y;
	mCv[3 * k + 2] = z;
}

void NurbsCurve3D::GenerateEv(double *ev, int point_count) const
{
	for (int k = 0; k < point_count; ++k) {
		double u = (double)k / (point_count - 1);
		POINT_3D p = NURBS_3D(mKnots, mKnotCount, mDegree, 0, u, (POINT_3D*)mCv, mCvCount);
		memcpy(ev + 3 * k, &p, sizeof(POINT_3D));
	}
}

NurbsCurve3D NurbsCurve3D::FromEv(const double *ev, int point_count, int degree)
{
    NurbsCurve3D inst(point_count, degree);
    double *uk_vector = new double[point_count];
    global_do_knots_chord3d((POINT_3D*)ev, point_count, uk_vector);
    global_do_knots(inst.GetKnots(), inst.GetKnotCount(), uk_vector, point_count, degree);

    int ret = nurbs_calc_control3d((POINT_3D*)ev, point_count, degree, uk_vector, inst.GetKnots(), inst.GetKnotCount(),
                                   (POINT_3D*)inst.GetCv());

    delete[]uk_vector;
    if (ret != E_OK) {
        //
        return NurbsCurve3D(); // Invalid
    }
    //
    return inst;
}
NurbsCurve3D *NurbsCurve3D::CreateFromEv(const double *ev, int point_count, int degree)
{
	NurbsCurve3D *inst = new NurbsCurve3D(point_count, degree);
	double *uk_vector = new double[point_count];
	global_do_knots_chord3d((POINT_3D*)ev, point_count, uk_vector);
	global_do_knots(inst->GetKnots(), inst->GetKnotCount(), uk_vector, point_count, degree);

	int ret = nurbs_calc_control3d((POINT_3D*)ev, point_count, degree, uk_vector, inst->GetKnots(), inst->GetKnotCount(),
		(POINT_3D*)inst->GetCv());

	delete[]uk_vector;
	if (ret != E_OK) {
		//
		delete inst;
		// No solution.
		return 0;
	}
	//
	return inst;
}
