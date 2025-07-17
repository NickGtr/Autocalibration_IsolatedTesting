#pragma once

#include <Eigen/QR>
#include <Eigen/Eigenvalues>

namespace rootba_povar {

    using Mat3 = Eigen::Matrix<double, 3, 3>;
    using Mat4 = Eigen::Matrix<double, 4, 4>;
    using Mat34 = Eigen::Matrix<double, 3, 4>;
    using Vec3 = Eigen::Matrix<double, 3, 1>;

    template <typename Scalar>
    void KRt_From_P(const Eigen::Matrix<Scalar, 3, 4> &P, Eigen::Matrix<Scalar, 3, 3> *K, Eigen::Matrix<Scalar, 3, 3> *R, Eigen::Matrix<Scalar, 3, 1> *t);

    template <typename Scalar>
    void K_From_ImageOfTheDualAbsoluteQuadric(const Eigen::Matrix<Scalar, 3, 3> &KKt, Eigen::Matrix<Scalar, 3, 3> *K);

    template <typename Scalar>
    void K_From_AbsoluteConic(const Eigen::Matrix<Scalar, 3, 3> &W, Eigen::Matrix<Scalar, 3, 3> *K);

    template <typename Scalar>
    class AutoCalibrationLinear {
    public:
        /** \brief Add a projection to be used for autocalibration.
         *
         *  \param P The projection matrix.
         *  \param width  The width of the image plane.
         *  \param height The height of the image plane.
         *
         *  The width and height parameters are used to normalize the projection
         *  matrix for improving numerical stability.  The don't need to be exact.
         */
        int AddProjection(const Eigen::Matrix<Scalar, 3, 4> &P, Scalar width, Scalar height);

        int AddProjection_1997paper(const Eigen::Matrix<Scalar, 3, 4> &P, Scalar width, Scalar height);

        /** \brief Computes the metric updating transformation.
         *
         *  \return The homography, H, that transforms the space into a metric space.
         *          If {P, X} is a projective reconstruction, then {P H, H^{-1} X} is
         *          a metric reconstruction.  Note that this follows the notation of
         *          HZ section 19.1 page 459, and not the notation of Pollefeys'
         *          paper [1].
         */
        Eigen::Matrix<Scalar, 4, 4> MetricTransformation(Eigen::Matrix<Scalar, 4, 4> *Q_final = nullptr, int *rank = nullptr);

    private:
        /** \brief Add constraints on the absolute quadric based assumptions on the
         *         parameters of one camera.
         *
         *  \param P The projection matrix of the camera in projective coordinates.
         */

        /** \brief Computes the constraint associated to elements of the DIAC.
         *
         *  \param P The projection used to project the absolute quadric.
         *  \param i Row of the DIAC.
         *  \param j Column of the DIAC.
         *  \return The coeficients of the element i, j of the dual image of the
         *          absolute conic when written as a linear combination of the
         *          elements of the absolute quadric.  There are 10 coeficients since
         *          the absolute quadric is represented by 10 numbers.
         */

        void AddProjectionConstraints(const Eigen::Matrix<Scalar, 3, 4> &P);
        void AddProjectionConstraints_1997paper(const Eigen::Matrix<Scalar, 3, 4> &P);

        Scalar Nullspace(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *A, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *nullspace);

        void SortEigenVectors(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &values,
                                     const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &vectors,
                                     Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *sorted_values,
                                     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *sorted_vectors);

        static Eigen::Matrix<Scalar, Eigen::Dynamic, 1> wc(const Eigen::Matrix<Scalar, 3, 4> &P, int i, int j);

        static Eigen::Matrix<Scalar, 4, 4> AbsoluteQuadricMatFromVec(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &q);

        static void NormalizeProjection(const Eigen::Matrix<Scalar, 3, 4> &P,
                                        Scalar width,
                                        Scalar height,
                                        Eigen::Matrix<Scalar, 3, 4> *P_new);

        static void DenormalizeProjection(const Eigen::Matrix<Scalar, 3, 4> &P,
                                          Scalar width,
                                          Scalar height,
                                          Eigen::Matrix<Scalar, 3, 4> *P_new);

        static void get_rid_of_principal_point(const Eigen::Matrix<Scalar, 3, 4> &P,
                                          Scalar width,
                                          Scalar height,
                                          Eigen::Matrix<Scalar, 3, 4> *P_new);

        static void get_back_principal_point(const Eigen::Matrix<Scalar, 3, 4> &P,
                                          Scalar width,
                                          Scalar height,
                                          Eigen::Matrix<Scalar, 3, 4> *P_new);

    private:
        std::vector<Eigen::Matrix<Scalar, 3, 4>> projections_; // The *normalized* projection matrices.
        std::vector<Scalar> widths_;
        std::vector<Scalar> heights_;
        std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> constraints_;  // Linear constraints on q.
    };
}
