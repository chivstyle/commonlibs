///
/// (c) 2019 chiv
///
#pragma once

namespace geometry {

    class Geometry {
    public:
        //
        Geometry(size_t number_of_vertices, size_t number_of_indices);
        Geometry(const Geometry& geometry);
        Geometry(Geometry&& geometry);
        Geometry();
        virtual ~Geometry();
        //
        struct Vertex {
            double x, y, z;
            //
            Vertex(double x_, double y_, double z_)
                : x(x_), y(y_), z(z_)
            {}
            Vertex() {}
        };
        //
        Geometry& operator=(const Geometry& geometry);
        Geometry& operator=(Geometry&& geometry);
        //
        bool IsValid() const;
        bool IsValidReally() const;
        //
        void Push(const Vertex& vertex);
        void Push(int vertex_index);
        void RemoveVertex(size_t idx = -1);
        void RemoveIndex(size_t idx = -1);
        // Resize capacity
        void ResizeIndices(size_t sz);
        void ResizeVertices(size_t sz);
        // SetVertex/SetIndex does not allocate memory, so before use them, please
        // allocate memory with resize, or create instance with "size" parameters.
        void SetVertex(size_t idx, const Vertex& vertex);
        void SetIndex(size_t idx, int vertex_index);
        // return vector.data()
        const int* GetIndices() const;
        // return vector.data()
        const Vertex* GetVertices() const;
        // return 0 if resource is empty
        size_t GetNumberOfIndices() const;
        // return 0 if resource is empty
        size_t GetNumberOfVertices() const;
        //
        bool SaveToSTL(const char* filepath);
    protected:
        class Resource;
        Resource* mRc;
    };


}
