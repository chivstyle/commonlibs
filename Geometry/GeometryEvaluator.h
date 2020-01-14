///
/// (c) 9102 chiv
///
#pragma once

namespace geometry {
    class Geometry;
    class GeometryEvaluator {
    public:
        GeometryEvaluator();
        virtual ~GeometryEvaluator();
        //
        struct GeometryMetrics {
            double SurfaceArea;
            double Volume;
        };
        GeometryMetrics Calculate(const Geometry& geometry);
    };


}

