///
/// (c) 2019 chiv
///
#include "Geometry.h"
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdint>

namespace geometry {

    class Geometry::Resource {
    public:
        Resource(size_t nv, size_t ni)
            : mIndices(ni)
            , mVertices(nv)
        {
        }
        Resource() {}
        //
        std::vector<int>& Indices() { return mIndices; }
        std::vector<Vertex>& Vertices() { return mVertices; }
        //
        void ResizeIndices(size_t sz)
        {
            mIndices.resize(sz);
        }
        void ResizeVertices(size_t sz)
        {
            mVertices.resize(sz);
        }
        //
    protected:
        std::vector<int> mIndices;
        std::vector<Vertex> mVertices;
    };

    Geometry::Geometry(size_t number_of_vertices, size_t number_of_indices)
        : mRc(new Resource(number_of_vertices, number_of_indices))
    {
    }
    Geometry::Geometry(const Geometry& geometry)
        : mRc(new Resource())
    {
        // copy.
        *mRc = *geometry.mRc;
    }
    Geometry::Geometry(Geometry&& geometry)
        : mRc(new Resource())
    {
        mRc->Indices() = std::move(geometry.mRc->Indices());
        mRc->Vertices() = std::move(geometry.mRc->Vertices());
    }
    Geometry::Geometry()
        : mRc(new Resource())
    {
    }
    Geometry::~Geometry()
    {
        delete mRc;
    }
    //
    void Geometry::ResizeIndices(size_t sz)
    {
        mRc->ResizeIndices(sz);
    }
    void Geometry::ResizeVertices(size_t sz)
    {
        mRc->ResizeVertices(sz);
    }
    //
    void Geometry::SetVertex(size_t idx, const Vertex& vertex)
    {
        assert(mRc != nullptr);
        mRc->Vertices()[idx] = vertex;
    }
    void Geometry::SetIndex(size_t idx, int vertex_index)
    {
        assert(mRc != nullptr);
        mRc->Indices()[idx] = vertex_index;
    }
    //
    Geometry& Geometry::operator=(const Geometry& geometry)
    {
        if (this != &geometry) {
            mRc->Indices() = geometry.mRc->Indices();
            mRc->Vertices() = geometry.mRc->Vertices();
        }
        return *this;
    }
    Geometry& Geometry::operator=(Geometry&& geometry)
    {
        if (this != &geometry) {
            mRc->Indices() = std::move(geometry.mRc->Indices());
            mRc->Vertices() = std::move(geometry.mRc->Vertices());
        }
        return *this;
    }
    //
    bool Geometry::IsValid() const
    {
        return !mRc->Indices().empty() && !mRc->Vertices().empty();
    }
    bool Geometry::IsValidReally() const
    {
        if (IsValid()) {
            for (size_t k = 0; k < mRc->Indices().size(); ++k) {
                if ((size_t)mRc->Indices()[k] >= mRc->Vertices().size()) return false;
            }
            return true;
        }
        return false;
    }
    //
    void Geometry::Push(const Vertex& vertex)
    {
        mRc->Vertices().push_back(vertex);
    }
    void Geometry::Push(int vertex_index)
    {
        mRc->Indices().push_back(vertex_index);
    }
    void Geometry::RemoveVertex(size_t idx)
    {
        if (idx == -1) {
            // remove last.
            mRc->Vertices().pop_back();
        }
        else {
            if (idx < mRc->Vertices().size()) {
                mRc->Vertices().erase(mRc->Vertices().begin() + idx);
            }
        }
    }
    void Geometry::RemoveIndex(size_t idx)
    {
        if (idx == -1) {
            // remove last.
            mRc->Indices().pop_back();
        } else {
            if (idx < mRc->Indices().size()) {
                mRc->Indices().erase(mRc->Indices().begin() + idx);
            }
        }
    }
    //
    const int* Geometry::GetIndices() const { return mRc->Indices().data(); }
    const Geometry::Vertex* Geometry::GetVertices() const { return mRc->Vertices().data(); }
    size_t Geometry::GetNumberOfIndices() const { return mRc->Indices().size(); }
    size_t Geometry::GetNumberOfVertices() const { return mRc->Vertices().size(); }
    //
    Geometry::Vertex operator-(const Geometry::Vertex& v1, const Geometry::Vertex& v2)
    {
        return Geometry::Vertex(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    //
    static __inline Geometry::Vertex _cross(const Geometry::Vertex& vec1, const Geometry::Vertex& vec2)
    {
        return Geometry::Vertex(
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x);
    }
    bool Geometry::SaveToSTL(const char* filepath)
    {
        FILE* fp = nullptr;
        errno_t er = fopen_s(&fp, filepath, "wb");
        if (er == 0) {
            char head[80];
            fwrite(head, 1, sizeof(head), fp);
            uint32_t triangles_count = (uint32_t)GetNumberOfIndices() / 3;
            fwrite(&triangles_count, 1, 4, fp);
            for (size_t k = 0; k < GetNumberOfIndices(); k += 3) {
                auto &v1 = GetVertices()[GetIndices()[k]];
                auto &v2 = GetVertices()[GetIndices()[k+1]];
                auto &v3 = GetVertices()[GetIndices()[k+2]];
#pragma pack(push)
#pragma pack(1)
                struct Node {
                    float normal[3];
                    float vert1[3];
                    float vert2[3];
                    float vert3[3];
                    short ab;
                } nb = {
                    { 0, 0, 1 },
                    { (float)v1.x, (float)v1.y, (float)v1.z },
                    { (float)v2.x, (float)v2.y, (float)v2.z },
                    { (float)v3.x, (float)v3.y, (float)v3.z },
                    0
                };
#if 0
                // generate normal
                auto normal = _cross(
                    Vertex(nb.vert2[0], nb.vert2[1], nb.vert2[2]) - Vertex(nb.vert1[0], nb.vert1[1], nb.vert1[2]),
                    Vertex(nb.vert3[0], nb.vert3[1], nb.vert3[2]) - Vertex(nb.vert1[0], nb.vert1[1], nb.vert1[2])
                    );
                memcpy(nb.normal, &normal, sizeof(normal));
#endif
                //
                fwrite(&nb, 1, sizeof(nb), fp);
            }
#pragma pack(pop)
            fclose(fp);
            //
            return true;
        }
        return false;
    }
}