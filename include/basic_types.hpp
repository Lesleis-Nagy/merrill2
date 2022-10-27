#pragma once

#include <array>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef WITH_MULTIPRECISION

#include "mpreal.h"

#endif // WITH_MULTIPRECISION

namespace merrill2 {

    /**
     * An implementation of a three dimensional cartesian vector.
     * @tparam Real the underlying data type for the calculation - usually ‘double’ or ‘mpreal’.
     */
    template<typename Real>
    class Vector3D {

    public:

        /**
         * Set the regularisation-epsilon value for **all** Vector3D objects of this type.
         * @param new_eps the new regularisation-epsilon.
         */
        static void set_eps(Real new_eps) {

            _eps = new_eps;
            _eps_squared = new_eps * new_eps;

        }

        /**
         * Retrieve the regularisation-epsilon.
         * @return the regularisation-epsilon.
         */
        static Real eps() {

            return _eps;

        }

        /**
         * Retrieve the regularisation-epsilon squared.
         * @return the regularisation-epsilon squared.
         */
        static Real eps_squared() {

            return _eps_squared;

        }

        /**
         * Create a three dimensional zero-vector object.
         */
        Vector3D() : _x(0), _y(0), _z(0) {}

        /**
         * Create a three dimensional vector object with the given x, y & z components along
         * with a regularisation-epsilon value.
         * @param x the vector x component.
         * @param y the vector y component.
         * @param z the vector z component.
         * @param eps the regularization-epsilon value.
         */
        Vector3D(Real x, Real y, Real z) : _x(std::move(x)), _y(std::move(y)), _z(std::move(z)) {}

        /**
         * Retrieve the vector's x-component.
         * @return the vector's x-component.
         */
        [[nodiscard]] inline Real x() const { return _x; }

        /**
         * Retrieve the vector's y-component.
         * @return the vector's y-component.
         */
        [[nodiscard]] inline Real y() const { return _y; }

        /**
         * Retrieve the vector's z-component.
         * @return the vector's z-component.
         */
        [[nodiscard]] inline Real z() const { return _z; }

    private:

        Real _x;
        Real _y;
        Real _z;

        [[maybe_unused]] static Real _eps;
        [[maybe_unused]] static Real _eps_squared;

    };

    // Initialize static eps & eps_squared values to defaults for double precision arithmetic.
    template<typename Real>
    Real Vector3D<Real>::_eps = 1E-7;

    template<typename Real>
    Real Vector3D<Real>::_eps_squared = 1E-14;

    /**
     * Redirection operator to display the vector.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param out the output stream.
     * @param v the vector to display.
     * @return the output stream with a representation of the input vector.
     */
    template<typename Real>
    std::ostream &operator<<(std::ostream &out, const Vector3D<Real> v) {

        out << "<" << v.x() << ", " << v.y() << ", " << v.z() << ">";
        return out;

    }

    /**
     * Vector addition operator.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param u the vector on the left hand side of the sum.
     * @param v the vector on the right hand side of the sum.
     * @return the sum of the two input vectors.
     */
    template<typename Real>
    Vector3D<Real> operator+(const Vector3D<Real> &u, const Vector3D<Real> &v) {

        return {u.x() + v.x(), u.y() + v.y(), u.z() + v.z()};

    }

    /**
     * Vector subtraction operator.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param u the vector on the left hand side of the operator.
     * @param v the vector on the right hand side of the operator.
     * @return the difference of two input vectors.
     */
    template<typename Real>
    Vector3D<Real> operator-(const Vector3D<Real> &u, const Vector3D<Real> &v) {

        return {u.x() - v.x(), u.y() - v.y(), u.z() - v.z()};

    }

    /**
     * Vector-scalar product operator.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param v the vector on the left hand side of the product.
     * @param lambda the scalar on the right hand side of the product.
     * @return the vector-scalar product.
     */
    template<typename Real>
    Vector3D<Real> operator*(const Vector3D<Real> &v, Real lambda) {

        return {lambda * v.x(), lambda * v.y(), lambda * v.z()};

    }

    /**
     * Scalar-vector product operator.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param lambda the scalar on the left hand side of the product.
     * @param v the vector on the right hand side of the product.
     * @return the scalar-vector product.
     */
    template<typename Real>
    Vector3D<Real> operator*(Real lambda, const Vector3D<Real> &v) {

        return {lambda * v.x(), lambda * v.y(), lambda * v.z()};

    }

    /**
     * Vector-scalar division.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param v the vector on the left hand side of the division.
     * @param lambda the scalar on the right hand side of the division.
     * @return the vector-scalar division.
     */
    template<typename Real>
    Vector3D<Real> operator/(const Vector3D<Real> &v, Real lambda) {

        return {v.x() / lambda, v.y() / lambda, v.z() / lambda};

    }

    /**
     * Vector dot product.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param u the vector on the left hand side of the dot product.
     * @param v the scalar on the right hand side of the dot product.
     * @return the vector dot product.
     */
    template<typename Real>
    Real dot(const Vector3D<Real> &u, const Vector3D<Real> &v) {

        return u.x() * v.x() + u.y() * v.y() + u.z() * v.z();

    }

    /**
     * Vector cross product.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param u the vector on the left hand side of the cross product.
     * @param v the scalar on the right hand side of the cross product.
     * @return the vector cross product.
     */template<typename Real>
    Vector3D<Real> cross(const Vector3D<Real> &u, const Vector3D<Real> &v) {

        return {u.y() * v.z() - u.z() * v.y(),
                -u.x() * v.z() + u.z() * v.x(),
                u.x() * v.y() - u.y() * v.x()};

    }

    /**
     * The norm of a vector
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param v the vector for which we seek the norm.
     * @return the norm of the input vector.
     */
    template<typename Real>
    Real norm(const Vector3D<Real> &v) {

        return sqrt(dot(v, v) + Vector3D<Real>::eps_squared());

    }

    /**
     * The norm squared of a vector.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param v the vector for which we seek the norm-squared.
     * @return the norm-squared of the input vector.
     */
    template<typename Real>
    Real norm_squared(const Vector3D<Real> &v) {

        return dot(v, v);

    }

    /**
     * Produce a normalised version of the input vector.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param v the vector that we seek to normalise.
     * @return the normalised input vector.
     */
    template<typename Real>
    Vector3D<Real> normalised(const Vector3D<Real> &v) {

        Real l = norm(v);
        return v / l;

    }

    ////////////////////////////////////// Geometry ////////////////////////////////////////////////////////////////////

    /**
     * Return the edge_length between two vector endpoints.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param lhs vector representing the start point of the edge.
     * @param rhs vector representing the end point of the edge.
     * @return the length of the edge.
     */
    template<typename Real>
    Real edge_length(const Vector3D<Real> &lhs, const Vector3D<Real> &rhs) {

        return norm(lhs - rhs);

    }

    /**
     * Return the edge center between two vector endpoints.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param r1 vector representing the start point of the edge.
     * @param r2 vector representing the end point of the edge.
     * @return the vector representing the center point of the edge.
     */
    template<typename Real>
    Vector3D<Real> edge_center(const Vector3D<Real> &r1, const Vector3D<Real> &r2) {

        return (r1 + r2) / 2.0;

    }

#ifdef WITH_MULTIPRECISION
    /**
     * Return the edge center between two vector endpoints - mpreal specific version.
     * @param r1 vector representing the start point of the edge.
     * @param r2 vector representing the end point of the edge.
     * @return the vector representing the center point of the edge.
     */
    template<>
    Vector3D<mpfr::mpreal> edge_center(const Vector3D<mpfr::mpreal> &r1, const Vector3D<mpfr::mpreal> &r2) {

        using mpfr::mpreal;

        return (r1 + r2) / mpreal(2.0);

    }
#endif // WITH_MULTIPRECISION

    /**
     * Return the orientation vector between two vector end points; this is the unit vector pointing from \f$r_1\f$ to
     * \f$r_2\f$.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param r1 vector representing the start point of the edge.
     * @param r2 vector representing the end point of the edge.
     * @return the unit vector pointing from \f$r_1\f$ to \f$r_2\f$.
     */
    template<typename Real>
    Vector3D<Real> edge_orientation(const Vector3D<Real> &r1, const Vector3D<Real> &r2) {

        return normalised(r2 - r1);

    }

    /**
     * Return the triangle normal vector assuming vertex clockwise winding \f$ r_1 \rightarrow r_2 \f$,
     * \f$ r_2 \rightarrow r_3 \f$ and \f$ r_3 \rightarrow r_1 \f$.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param r1 vector representing a point on the triangle.
     * @param r2 vector representing a point on the triangle.
     * @param r3 vector representing a point on the triangle.
     * @return the triangle normal vector.
     */
    template<typename Real>
    Vector3D<Real> triangle_normal(const Vector3D<Real> &r1, const Vector3D<Real> &r2, const Vector3D<Real> &r3) {

        return normalised(cross(r2 - r1, r3 - r1));

    }

    /**
     * Return the triangle center vector.
     * @tparam Real the underlying data type for the calculation - usually 'double' or 'mpreal'.
     * @param r1 vector representing a point on the triangle.
     * @param r2 vector representing a point on the triangle.
     * @param r3 vector representing a point on the triangle.
     * @return the triangle center vector.
     */
    template<typename Real>
    Vector3D<Real> triangle_center(const Vector3D<Real> &r1, const Vector3D<Real> &r2, const Vector3D<Real> &r3) {

        Vector3D<Real> sum = (r1 + r2) + r3;
        return sum / 3.0;

    }

#ifdef WITH_MULTIPRECISION
    /**
     * Return the triangle center vector - mpreal specific version.
     * @param r1 vector representing a point on the triangle.
     * @param r2 vector representing a point on the triangle.
     * @param r3 vector representing a point on the triangle.
     * @return the triangle center vector.
     */
    template<>
    Vector3D<mpfr::mpreal> triangle_center(const Vector3D<mpfr::mpreal> &r1,
                                           const Vector3D<mpfr::mpreal> &r2,
                                           const Vector3D<mpfr::mpreal> &r3) {

        using mpfr::mpreal;

        Vector3D<mpreal> sum = (r1 + r2) + r3;
        return sum / mpreal(3.0);

    }
#endif // WITH_MULTIPRECISION

    template<typename Real>
    Vector3D<Real> tetrahedron_center(const Vector3D<Real> &r1,
                                      const Vector3D<Real> &r2,
                                      const Vector3D<Real> &r3,
                                      const Vector3D<Real> &r4) {

        Vector3D<Real> sum = ((r1 + r2) + r3) + r4;
        return sum / 4.0;

    }

#ifdef WITH_MULTIPRECISION
    /**
     * Return the tetrahedron center vector - mpreal specific version.
     * @param r1 vector representing a point on the tetrahedron.
     * @param r2 vector representing a point on the tetrahedron.
     * @param r3 vector representing a point on the tetrahedron.
     * @param r4 vector representing a point on the tetrahedron.
     * @return the tetrahedron center vector.
     */
    template<>
    Vector3D<mpfr::mpreal> tetrahedron_center(const Vector3D<mpfr::mpreal> &r1,
                                              const Vector3D<mpfr::mpreal> &r2,
                                              const Vector3D<mpfr::mpreal> &r3,
                                              const Vector3D<mpfr::mpreal> &r4) {

        using mpfr::mpreal;

        Vector3D<mpreal> sum = ((r1 + r2) + r3) + r4;
        return sum / mpreal(4.0);

    }
#endif // WITH_MULTIPRECISION

    template<typename Real>
    using Vector3dList = std::vector<Vector3D<Real>>;

    template<typename Integral>
    class PointIndex {

    public:
        explicit PointIndex(Integral n0) : _n0{n0} {}

        Integral operator()() const { return _n0; }

        void add_subdomain_index(Integral i) { _subdomain_indices.insert(i); }

        const std::set<Integral> &subdomain_indices() const { return _subdomain_indices; }

    private:

        Integral _n0;

        std::set<Integral> _subdomain_indices;

    };

    template<typename Integral>
    class LineIndex {

    public:
        LineIndex(Integral n0, Integral n1) : _n{n0, n1} {}

        Integral operator()(size_t i) const { return _n[i % 2]; }

        PointIndex<Integral> contract(size_t i) {
            switch (i) {
                case 0:
                    return PointIndex{_n[1]};
                case 1:
                    return PointIndex{_n[0]};
                default:
                    throw std::runtime_error("Cannot perform contraction on LineIndex - out of bounds.");
            }
        }

        void add_subdomain_index(Integral i) { _subdomain_indices.insert(i); }

        const std::set<Integral> &subdomain_indices() const { return _subdomain_indices; }

    private:

        std::array<Integral, 2> _n;

        std::set<Integral> _subdomain_indices;

    };

    template<typename Integral>
    class TriangleIndex {

    public:

        TriangleIndex(Integral n0, Integral n1, Integral n2) : _n{n0, n1, n2} {}

        Integral operator()(size_t i) const { return _n[i % 3]; }

        LineIndex<Integral> con(size_t i) {
            switch (i) {
                case 0:
                    return LineIndex{_n[1], _n[2]};
                case 1:
                    return LineIndex{_n[0], _n[2]};
                case 2:
                    return LineIndex{_n[0], _n[1]};
                default:
                    throw std::runtime_error("Cannot perform contraction on TriangleIndex - out of bounds.");
            }
        }

        template<typename Real>
        Vector3D<Real> normal(const std::vector<Vector3D<Real>> &verts) const {

            auto v0 = verts[0];
            auto v1 = verts[1];
            auto v2 = verts[2];

            return normalized(cross(v0 - v1, v2 - v1));

        }

        void add_subdomain_index(Integral i) { _subdomain_indices.insert(i); }

        const std::set<Integral> &subdomain_indices() const { return _subdomain_indices; }

    private:

        std::array<Integral, 3> _n;

        std::set<Integral> _subdomain_indices;

    };

    template<typename Integral>
    class TetrahedralIndex {

    public:

        template<typename Real>
        static TetrahedralIndex
        create(Integral n0, Integral n1, Integral n2, Integral n3, Integral subdomain_index,
               const Vector3dList<Real> &vertices) {

            auto v0 = vertices[n0];
            auto v1 = vertices[n1];
            auto v2 = vertices[n2];
            auto v3 = vertices[n3];

            // Calculate triangle center for some triangle (say the con(3) triangle).
            auto tr_c = triangle_center(v0, v1, v2);

            // Calculate the tetrahedron center.
            auto te_c = tetrahedron_center(v0, v1, v2, v3);

            // Calculate the normalized center to center vector.
            auto c2c = normalised(tr_c - te_c);

            // Calculate the test triangle normal.
            auto n = normalized(cross(v0 - v1, v2 - v1));

            // Make index correction depending on the dot product of `c2x` and `n` vectors.
            if (dot(c2c, n) < 0) {
                // The `c2c` and `n` vectors point in opposite directions to each other (w.r.t. the triangle <v0, v1, v2>),
                // therefore the tetrahedron has incorrect winding.
                return TetrahedralIndex(n1, n0, n2, n3, subdomain_index);
            } else {
                // The `c2c` and `n` vectors point in similar directions to each other (w.r.t. the triangle <v0, v1, v2>),
                // therefore the tetrahedron has correct winding.
                return TetrahedralIndex(n0, n1, n2, n3, subdomain_index);
            }

        }

        Integral operator[](size_t i) { return _n[i % 4]; }

        TriangleIndex<Integral> con(size_t i) {
            switch (i) {
                case 0:
                    return TriangleIndex{_n[1], _n[3], _n[2]};
                case 1:
                    return TriangleIndex{_n[0], _n[2], _n[3]};
                case 2:
                    return TriangleIndex{_n[0], _n[3], _n[1]};
                case 3:
                    return TriangleIndex{_n[0], _n[1], _n[2]};
                default:
                    throw std::runtime_error("Cannot perform contraction on TriangleIndex - out of bounds.");
            }
        }

        void set_subdomain_index(Integral i) { _subdomain_index = i; }

        Integral subdomain_index() const { return _subdomain_index; }

    private:

        TetrahedralIndex(Integral n0, Integral n1, Integral n2, Integral n3, Integral subdomain_index) :
                _n{n0, n1, n2, n3},
                _subdomain_index(subdomain_index) {}

        std::array<Integral, 4> _n;

        Integral _subdomain_index;

    };

} // merrill2