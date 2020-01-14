///
/// (c) 2019 chiv
///
#include "GeometryEvaluator.h"
#include "Geometry.h"
#include <cmath>
///
namespace geometry {

    GeometryEvaluator::GeometryEvaluator() {}
    GeometryEvaluator::~GeometryEvaluator() {}
    //
    static __inline double _volume(const Geometry::Vertex& p1, const Geometry::Vertex&  p2, const Geometry::Vertex&  p3)
    {
        double v321 = p3.x * p2.y * p1.z;
        double v231 = p2.x * p3.y * p1.z;
        double v312 = p3.x * p1.y * p2.z;
        double v132 = p1.x * p3.y * p2.z;
        double v213 = p2.x * p1.y * p3.z;
        double v123 = p1.x * p2.y * p3.z;
        // v123 == v231 == v312 == -132 == -213 == -321
        return (-v321 + v231 + v312 - v132 - v213 + v123) / 6;
    }
    //
    GeometryEvaluator::GeometryMetrics GeometryEvaluator::Calculate(const Geometry& geometry)
    {
        double s = 0, v = 0;
        for (size_t k = 0; k < geometry.GetNumberOfIndices(); k += 3) {
            const Geometry::Vertex vertex1 = geometry.GetVertices()[geometry.GetIndices()[k]];
            const Geometry::Vertex vertex2 = geometry.GetVertices()[geometry.GetIndices()[k+1]];
            const Geometry::Vertex vertex3 = geometry.GetVertices()[geometry.GetIndices()[k+2]];
            double a = sqrt(pow(vertex1.x - vertex2.x, 2) + pow(vertex1.y - vertex2.y, 2) + pow(vertex1.z - vertex2.z, 2));
            double b = sqrt(pow(vertex2.x - vertex3.x, 2) + pow(vertex2.y - vertex3.y, 2) + pow(vertex2.z - vertex3.z, 2));
            double c = sqrt(pow(vertex3.x - vertex1.x, 2) + pow(vertex3.y - vertex1.y, 2) + pow(vertex3.z - vertex1.z, 2));
            double p = (a + b + c) / 2;
            s += sqrt(p * (p - a) * (p - b) * (p - c));
            //
            v += _volume(vertex1, vertex2, vertex3);
        }
        return GeometryMetrics{ s, v };
    }
}