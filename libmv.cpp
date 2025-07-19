//
// Created by Simon on 14/03/2025.
//

#include "libmv.hpp"
#include <iostream>

namespace rootba_povar {


    // Nicolas
    template <typename Scalar>
    void K_From_ImageOfTheDualAbsoluteQuadric(const Eigen::Matrix<Scalar, 3, 3> &KKt, Eigen::Matrix<Scalar, 3, 3> *K) {
        // Does the same thing as K_From_AbsoluteConic but is properly named and stays in the dual space : no need to dualize (=take an inverse).
        // Use an upper triangular Choleski to uncover K from KKt


        Eigen::Matrix<Scalar, 3, 3> flipped_KKt;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                flipped_KKt(i,j) = KKt(2 - i, 2 - j);
            }
        }

        Eigen::LLT<Eigen::Matrix<Scalar, 3, 3>> llt(flipped_KKt);
        Eigen::Matrix<Scalar, 3, 3> L = llt.matrixL();

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                (*K)(i,j) = L(2 - i, 2 - j);
            }
        }

        // Resolve sign ambiguities assuming positive diagonal.
        for (int j = 0; j < 3; ++j) {
            if ((*K)(j, j) < 0) {
                for (int i = 0; i < 3; ++i) {
                    (*K)(i, j) = -(*K)(i, j);
                }
            }
        }
    }


    //@Simon: from llvm library: https://github.com/libmv/libmv/blob/master/src/libmv/multiview/autocalibration.cc

    template <typename Scalar>
    void K_From_AbsoluteConic(const Eigen::Matrix<Scalar, 3, 3> &W, Eigen::Matrix<Scalar, 3, 3> *K) {
        // To compute upper-triangular Cholesky, we flip the indices of the input
        // matrix, compute lower-triangular Cholesky, and then unflip the result.
        // Nicolas : K is not recovered from the absolute conic here but from the projection of the absolute quadric

        Eigen::Matrix<Scalar, 3, 3> dual = W.inverse();
        Eigen::Matrix<Scalar, 3, 3> flipped_dual;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                flipped_dual(i,j) = dual(2 - i, 2 - j);
            }
        }

        Eigen::LLT<Eigen::Matrix<Scalar, 3, 3>> llt(flipped_dual);
        Eigen::Matrix<Scalar, 3, 3> L = llt.matrixL();

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                (*K)(i,j) = L(2 - i, 2 - j);
            }
        }

        // Resolve sign ambiguities assuming positive diagonal.
        for (int j = 0; j < 3; ++j) {
            if ((*K)(j, j) < 0) {
                for (int i = 0; i < 3; ++i) {
                    (*K)(i, j) = -(*K)(i, j);
                }
            }
        }
    }

    // llvm library, see https://github.com/libmv/libmv/blob/master/src/libmv/multiview/projection.cc
    template <typename Scalar>
    void KRt_From_P(const Eigen::Matrix<Scalar, 3, 4> &P, Eigen::Matrix<Scalar, 3, 3> *Kp, Eigen::Matrix<Scalar, 3, 3> *Rp, Eigen::Matrix<Scalar, 3, 1> *tp) {
        // Decompose using the RQ decomposition HZ A4.1.1 pag.579.
        Eigen::Matrix<Scalar, 3, 3> K = P.block(0, 0, 3, 3);

        Eigen::Matrix<Scalar, 3, 3> Q;
        Q.setIdentity();

        // Set K(2,1) to zero.
        if (K(2,1) != 0) {
            Scalar c = -K(2,2);
            Scalar s = K(2,1);
            Scalar l = sqrt(c * c + s * s);
            c /= l; s /= l;
            Eigen::Matrix<Scalar, 3, 3> Qx;
            Qx << 1, 0, 0,
                    0, c, -s,
                    0, s, c;
            K = K * Qx;
            Q = Qx.transpose() * Q;
        }
        // Set K(2,0) to zero.
        if (K(2,0) != 0) {
            Scalar c = K(2,2);
            Scalar s = K(2,0);
            Scalar l = sqrt(c * c + s * s);
            c /= l; s /= l;
            Eigen::Matrix<Scalar, 3, 3> Qy;
            Qy << c, 0, s,
                    0, 1, 0,
                    -s, 0, c;
            K = K * Qy;
            Q = Qy.transpose() * Q;
        }
        // Set K(1,0) to zero.
        if (K(1,0) != 0) {
            Scalar c = -K(1,1);
            Scalar s = K(1,0);
            Scalar l = sqrt(c * c + s * s);
            c /= l; s /= l;
            Eigen::Matrix<Scalar, 3, 3> Qz;
            Qz << c,-s, 0,
                    s, c, 0,
                    0, 0, 1;
            K = K * Qz;
            Q = Qz.transpose() * Q;
        }

        Eigen::Matrix<Scalar, 3, 3> R = Q;

        // Ensure that the diagonal is positive.
        // TODO(pau) Change this to ensure that:
        //  - K(0,0) > 0
        //  - K(2,2) = 1
        //  - det(R) = 1
        if (K(2,2) < 0) {
            K = -K;
            R = -R;
        }
        if (K(1,1) < 0) {
            Eigen::Matrix<Scalar, 3, 3> S;
            S << 1, 0, 0,
                    0,-1, 0,
                    0, 0, 1;
            K = K * S;
            R = S * R;
        }
        if (K(0,0) < 0) {
            Eigen::Matrix<Scalar, 3, 3> S;
            S << -1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;
            K = K * S;
            R = S * R;
        }

        // Compute translation.
        Eigen::Matrix<Scalar, 3, 1> p(3);
        p << P(0,3), P(1,3), P(2,3);
        // TODO(pau) This sould be done by a SolveLinearSystem(A, b, &x) call.
        // TODO(keir) use the eigen LU solver syntax...
        Eigen::Matrix<Scalar, 3, 1> t = K.inverse() * p;

        // scale K so that K(2,2) = 1
        K = K / K(2,2);

        *Kp = K;
        *Rp = R;
        *tp = t;
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::NormalizeProjection(const Eigen::Matrix<Scalar, 3, 4> &P,
                             Scalar width,
                             Scalar height,
                             Eigen::Matrix<Scalar, 3, 4> *P_new) {
        Eigen::Matrix<Scalar, 3, 3> T;
        T << width + height,              0,  width / 2,
                0, width + height, height / 2,
                0,              0,                1;
        *P_new = T.inverse() * P;
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::get_rid_of_principal_point(const Eigen::Matrix<Scalar, 3, 4> &P,
                            Scalar width,
                            Scalar height,
                            Eigen::Matrix<Scalar, 3, 4> *P_new){
        // Nicolas : premultiply the projection matrix to get rid of the principal point assuming that it is dead center.
        Eigen::Matrix<Scalar, 3, 3> T;
        T << 1, 0, - width / 2,
             0, 1,  -height / 2,
             0, 0, 1;
        *P_new = T * P;                        
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::get_back_principal_point(const Eigen::Matrix<Scalar, 3, 4> &P,
                                                      Scalar width,
                                                      Scalar height,
                                                      Eigen::Matrix<Scalar, 3, 4> *P_new) {
        // Nicolas : postmultiply the projection matrix to recover the original principal point after the calculation
        Eigen::Matrix<Scalar, 3, 3> T;
        T << 1, 0, width / 2,
             0, 1, height / 2,
             0, 0, 1;
        *P_new = T * P;
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::DenormalizeProjection(const Eigen::Matrix<Scalar, 3, 4> &P,
                                                      Scalar width,
                                                      Scalar height,
                                                      Eigen::Matrix<Scalar, 3, 4> *P_new) {
        Eigen::Matrix<Scalar, 3, 3> T;
        T << width + height,              0,  width / 2,
                0, width + height, height / 2,
                0,              0,                1;
        *P_new = T * P;
    }

    template <typename Scalar>
    Eigen::Matrix<Scalar, 4, 4> AutoCalibrationLinear<Scalar>::AbsoluteQuadricMatFromVec(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &q) {
        Eigen::Matrix<Scalar, 4, 4> Q;
        Q << q(0), q(1), q(2), q(3),
                q(1), q(4), q(5), q(6),
                q(2), q(5), q(7), q(8),
                q(3), q(6), q(8), q(9);
        return Q;
    }

    template <typename Scalar>
    Scalar AutoCalibrationLinear<Scalar>::Nullspace(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *A, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *nullspace) {
        Eigen::JacobiSVD<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> svd(*A, Eigen::ComputeFullV);
        (*nullspace) = svd.matrixV().col(A->cols()-1);
        if (A->rows() >= A->cols())
            return svd.singularValues()(A->cols()-1);
        else
            return 0.0;
    }

    template <typename Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> AutoCalibrationLinear<Scalar>::wc(const Eigen::Matrix<Scalar, 3, 4> &P, int i, int j) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> constraint(10);
        for (int k = 0; k < 10; ++k) {
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> q = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(10);
            q(k) = 1;
            Eigen::Matrix<Scalar, 4, 4> Q = AbsoluteQuadricMatFromVec(q);

            Eigen::Matrix<Scalar, 3, 3> w = P * Q * P.transpose();

            constraint(k) = w(i, j);
        }
        return constraint;
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::AddProjectionConstraints(const Eigen::Matrix<Scalar, 3, 4> &P) {
        Scalar nu = 1.0;

        // Non-extreme focal lenght.
        constraints_.push_back((wc(P, 0, 0) - wc(P, 2, 2)) / 9 / nu);
        constraints_.push_back((wc(P, 1, 1) - wc(P, 2, 2)) / 9 / nu);

        // Aspect ratio is near 1.
        constraints_.push_back((wc(P, 0, 0) - wc(P, 1, 1)) / 0.2 / nu);

        // No skew and principal point near 0,0.
        // Note that there is a typo in the Pollefeys' paper: the 0.01 is not at the
        // correct equation.
        constraints_.push_back(wc(P, 0, 1) / 0.01 / nu);
        constraints_.push_back(wc(P, 0, 2) / 0.1 / nu);
        constraints_.push_back(wc(P, 1, 2) / 0.1 / nu);
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::AddProjectionConstraints_1997paper(const Eigen::Matrix<Scalar, 3, 4> &P) {
        Scalar nu = 1.0;

        // Aspect ratio is near 1.
        constraints_.push_back((wc(P, 0, 0) - wc(P, 1, 1)) / nu);

        // No skew and principal point near 0,0.
        // Note that there is a typo in the Pollefeys' paper: the 0.01 is not at the
        // correct equation.
        constraints_.push_back(wc(P, 0, 1) / nu);
        constraints_.push_back(wc(P, 0, 2) / nu);
        constraints_.push_back(wc(P, 1, 2) / nu);
    }

    template <typename Scalar>
    void AutoCalibrationLinear<Scalar>::SortEigenVectors(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &values,
                                 const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &vectors,
                                 Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *sorted_values,
                                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *sorted_vectors) {
        // Compute eigenvalues order.
        std::pair<Scalar, int> order[4];
        for (int i = 0; i < 4; ++i) {
            order[i].first = -values(i);
            order[i].second = i;
        }
        std::sort(order, order + 4);

        sorted_values->resize(4);
        sorted_vectors->resize(4,4);
        for (int i = 0; i < 4; ++i) {
            (*sorted_values)(i) = values[order[i].second];
            sorted_vectors->col(i) = vectors.col(order[i].second);
        }
    }

    template <typename Scalar>
    int AutoCalibrationLinear<Scalar>::AddProjection(const Eigen::Matrix<Scalar, 3, 4> &P,
                                                     Scalar width, Scalar height) {
        Eigen::Matrix<Scalar, 3, 4> P_normalized;
        //P_normalized = P;
        NormalizeProjection(P, width, height, &P_normalized);

        AddProjectionConstraints(P_normalized);

        // Store input
        projections_.push_back(P_normalized);
        widths_.push_back(width);
        heights_.push_back(height);

        return projections_.size() - 1;
    }

    template <typename Scalar>
    int AutoCalibrationLinear<Scalar>::AddProjection_1997paper(const Eigen::Matrix<Scalar, 3, 4> &P,
                                                     Scalar width, Scalar height) {
        Eigen::Matrix<Scalar, 3, 4> P_normalized;
        //P_normalized = P;
        get_rid_of_principal_point(P, width, height, &P_normalized);

        AddProjectionConstraints_1997paper(P_normalized);

        // Store input
        projections_.push_back(P_normalized);
        widths_.push_back(width);
        heights_.push_back(height);

        return projections_.size() - 1;
    }


    template <typename Scalar>
    Eigen::Matrix<Scalar, 4, 4> AutoCalibrationLinear<Scalar>::MetricTransformation(Eigen::Matrix<Scalar, 4, 4> *Q_final, int* rank) {
        // Compute the dual absolute quadric, Q.
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> A(constraints_.size(), 10);
        for (int i = 0; i < A.rows(); ++i) {
            A.row(i) = constraints_[i];
        }
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> q;
        Nullspace(&A, &q);
        Eigen::Matrix<Scalar, 4, 4> Q = AbsoluteQuadricMatFromVec(q);
        if (Q_final){
            *Q_final = Q;
        }
        // TODO(pau) force rank 3.

        // Compute a transformation to a metric frame by decomposing Q.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 4, 4> > eigen_solver(Q);

        // Eigen values should be possitive,
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> temp_values = eigen_solver.eigenvalues();
        if (temp_values.sum() < 0) {
            temp_values = -temp_values;
        }

        // and sorted, so that last one is 0.
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> eigenvalues;
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> eigenvectors;
        SortEigenVectors(temp_values, eigen_solver.eigenvectors(),
                         &eigenvalues, &eigenvectors);

        if (rank) {
            double relative_eps = 1e-3;
            double threshold = relative_eps * eigenvalues(0);
            *rank = 0;
            for (int i = 0; i < 4; ++i) {
                if (eigenvalues(i) > threshold) (*rank)++;
            }
            // if (*rank !=3){
            //     std::cout << "Rank : " << *rank << "Eigenvalues : " << eigenvalues.transpose() << std::endl;
            // }
        }

        // Compute the transformation from the eigen descomposition.  See last
        // paragraph of page 3 in
        //   "Autocalibration and the absolute quadric" by B. Triggs.
        eigenvalues(3) = 1;
        eigenvalues = eigenvalues.array().sqrt();
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> H = eigenvectors * eigenvalues.asDiagonal();
        return H;
    }
#ifdef ROOTBA_INSTANTIATIONS_FLOAT
    template void K_From_ImageOfTheDualAbsoluteQuadric(const Eigen::Matrix<float, 3, 3> &KKt, Eigen::Matrix<float, 3, 3> *K);
    template void K_From_AbsoluteConic(const Eigen::Matrix<float, 3, 3> &W, Eigen::Matrix<float, 3, 3> *K);
    template void KRt_From_P(const Eigen::Matrix<float, 3, 4> &P, Eigen::Matrix<float, 3, 3> *Kp, Eigen::Matrix<float, 3, 3> *Rp, Eigen::Matrix<float, 3, 1> *tp);
    template class AutoCalibrationLinear<float>;
#endif

#ifdef ROOTBA_INSTANTIATIONS_DOUBLE
    template void K_From_ImageOfTheDualAbsoluteQuadric(const Eigen::Matrix<double, 3, 3> &KKt, Eigen::Matrix<double, 3, 3> *K);
    template void K_From_AbsoluteConic(const Eigen::Matrix<double, 3, 3> &W, Eigen::Matrix<double, 3, 3> *K);
    template void KRt_From_P(const Eigen::Matrix<double, 3, 4> &P, Eigen::Matrix<double, 3, 3> *Kp, Eigen::Matrix<double, 3, 3> *Rp, Eigen::Matrix<double, 3, 1> *tp);
    template class AutoCalibrationLinear<double>;
#endif
}