#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>

using namespace boost;
using namespace boost::filesystem;

#include "casadi.hpp"
#include "orderparameter.hpp"

template SX energy<SX>(SX& fin, SX& J, SX& U0, SX& dU, SX mu);
template SX canonical<SX>(SX& fin, SX& J, SX& U0, SX& dU, SX mu);
template SX energy<double>(SX& fin, SX& J, SX& U0, SX& dU, double mu);
template SX canonical<double>(SX& fin, SX& J, SX& U0, SX& dU, double mu);

template SX energy2<SX>(SX& fin, SX& J, SX& U0, SX& dU, SX mu);
template SX energy2<double>(SX& fin, SX& J, SX& U0, SX& dU, double mu);

//double JW(double W) {
//    return alpha * (W * W) / (Ng * Ng + W * W);
//}
//
//double JWij(double Wi, double Wj) {
//    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
//}
//
//double UW(double W) {
//    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
//}
//
//SX JW(SX W) {
//    return alpha * (W * W) / (Ng * Ng + W * W);
//}
//
//SX JWij(SX Wi, SX Wj) {
//    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
//}
//
//SX UW(SX W) {
//    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
//}
//

complex<double> dot(vector<complex<double>>&v, vector<complex<double>>&w) {
    complex<double> res = 0;
    for (int i = 0; i < v.size(); i++) {
        res += ~v[i] * w[i];
    }
    return res;
}

namespace casadi {

    inline bool isnan(SX& sx) {
        return sx.at(0).isNan();
    }

    inline bool isinf(SX sx) {
        return sx.at(0).isInf();
    }
}

double DynamicsProblem::U00;
vector<double> DynamicsProblem::J0;

vector<double> DynamicsProblem::x0;
vector<vector<complex<double>>> DynamicsProblem::f0;

void DynamicsProblem::setup(double W, double mu, vector<double>& xi) {

    SX f = SX::sym("f", 2 * L * dim);
    SX dU = SX::sym("dU", L);
    SX J = SX::sym("J", L);
    SX U0 = SX::sym("U0");

    SX Wi = W;

    U0 = UW(W);
    for (int i = 0; i < L; i++) {
        J[i] = JWij(W * xi[i], W * xi[mod(i + 1)]);
        dU[i] = UW(W * xi[i]) - U0;
    }

    SX E = energy(f, J, U0, dU, mu);

    SXFunction nlp("nlp", nlpIn("x", f), nlpOut("f", E));
    NlpSolver solver("solver", "ipopt", nlp, make_dict("hessian_approximation", "limited-memory", "linear_solver", "ma86", "print_level", 0, "print_time", false));

    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> uni(-1, 1);

    vector<double> xrand(2 * L*dim, 1);
    rng.seed();
    for (int i = 0; i < 2 * L * dim; i++) {
        xrand[i] = uni(rng);
    }

    map<string, DMatrix> arg;
    arg["lbx"] = -1;
    arg["ubx"] = 1;
    arg["x0"] = xrand;

    map<string, DMatrix> res = solver(arg);
    x0 = res["x"].nonzeros();

    vector<complex<double>> f0i(dim);
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            f0i[n] = complex<double>(x0[2 * (i * dim + n)], x0[2 * (i * dim + n) + 1]);
        }
        double nrm = sqrt(abs(dot(f0i, f0i)));
        for (int n = 0; n <= nmax; n++) {
            x0[2 * (i * dim + n)] /= nrm;
            x0[2 * (i * dim + n) + 1] /= nrm;
        }
    }

    U00 = UW(W);
    J0 = vector<double>(L);
    for (int i = 0; i < L; i++) {
        J0[i] = JWij(W * xi[i], W * xi[mod(i + 1)]);
    }

    f0 = vector<vector<complex<double>>>(L, vector<complex<double>>(dim));
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            f0[i][n] = complex<double>(x0[2 * (i * dim + n)], x0[2 * (i * dim + n) + 1]);
        }
    }

}

DynamicsProblem::DynamicsProblem(double Wi, double Wf, double mu, vector<double> xi) {

    SX f = SX::sym("f", 2 * L * dim);
    SX dU = SX::sym("dU", L);
    SX J = SX::sym("J", L);
    SX U0 = SX::sym("U0");
    SX t = SX::sym("t");

    SX tau = SX::sym("tau");

    SX Wt = if_else(t < tau, Wi + (Wf - Wi) * t / tau, Wf + (Wi - Wf) * (t - tau) / tau);

    U0 = UW(Wt);
    for (int i = 0; i < L; i++) {
        J[i] = JWij(Wt * xi[i], Wt * xi[mod(i + 1)]);
        dU[i] = UW(Wt * xi[i]) - U0;
    }

    SX E = energy(f, J, U0, dU, mu);
    SXFunction Ef = SXFunction("E", {f, t, tau} , {E});
    E0 = SXFunction("E0",{f}, Ef(vector<SX>{f, 0, 1}));

    SX S = canonical(f, J, U0, dU, mu);
    SXFunction St("St", {t}, {S});
    SX Sdt = St.gradient()(vector<SX>{t})[0];


    SXFunction HSr("HSr", {f}, {Sdt});
    SX HSrdf = HSr.gradient()(vector<SX>{f})[0];
    SXFunction HSi("HSi", {f}, {-E});
    SX HSidf = HSi.gradient()(vector<SX>{f})[0];

    SX ode = SX::sym("ode", 2 * L * dim);
    for (int j = 0; j < L * dim; j++) {
        ode[2 * j] = 0;
        ode[2 * j + 1] = 0;
        try {
            ode[2 * j] += 0.5 * HSrdf[2 * j];
        }
        catch (CasadiException& e) {
        }
        try {
            ode[2 * j] -= 0.5 * HSidf[2 * j + 1];
        }
        catch (CasadiException& e) {
        }
        try {
            ode[2 * j + 1] += 0.5 * HSidf[2 * j];
        }
        catch (CasadiException& e) {
        }
        try {
            ode[2 * j + 1] += 0.5 * HSrdf[2 * j + 1];
        }
        catch (CasadiException& e) {
        }
    }
    ode_func = SXFunction("ode", daeIn("t", t, "x", f, "p", tau), daeOut("ode", ode));

//    double tauf = 2e-7;
//    double dt = 0.9e-9;
//    //    integrator = Integrator("integrator", "rk", ode_func, make_dict("t0", 0, "tf", 2*tauf, "number_of_finite_elements", ceil((2*tauf)/dt)));
//    integrator = Integrator("integrator", "cvodes", ode_func, make_dict("t0", 0, "tf", 3 * tauf, "exact_jacobian", false, "max_num_steps", 100000));
    Integrator integrator("integrator", "cvodes", ode_func, make_dict("exact_jacobian", false));
}

//void DynamicsProblem::evolve(double tau, results& result) {
//    double tauf = 2e-6;
//    double dt = 0.9e-9;
////    Integrator integrator("integrator", "rk", ode_func, make_dict("t0", 0, "tf", 2*tau, "number_of_finite_elements", ceil((2*tauf)/dt)));
//    Integrator integrator("integrator", "cvodes", ode_func, make_dict("t0", 0, "tf", 2 * tau, "exact_jacobian", false, "max_num_steps", 100000));
//
//    ptime start_time = microsec_clock::local_time();
//
//    result.f0 = f0;
//
//    map<string, DMatrix> res = integrator(make_map("x0", DMatrix(x0), "p", {tau}));
//    vector<double> xf = res["xf"].nonzeros();
//
//    ff = vector<vector<complex<double>>>(L, vector<complex<double>>(dim));
//    for (int i = 0; i < L; i++) {
//        for (int n = 0; n <= nmax; n++) {
//            ff[i][n] = complex<double>(xf[2 * (i * dim + n)], xf[2 * (i * dim + n) + 1]);
//        }
//        double nrm = sqrt(abs(dot(ff[i], ff[i])));
//        for (int n = 0; n <= nmax; n++) {
//            ff[i][n] /= nrm;
//        }
//    }
//    result.ff = ff;
//
//    result.J0 = J0;
//    result.U0 = U00;
//
//    result.b0 = vector<complex<double>>(L);
//    result.bf = vector<complex<double>>(L);
//    for (int i = 0; i < L; i++) {
//        result.b0[i] = b(f0, i, J0, U00);
//        result.bf[i] = b(ff, i, J0, U00);
//    }
//
//    vector<double> grad;
//    result.Ei = E0(vector<DMatrix>{x0})[0].toScalar();
//    result.Ef = E0(vector<DMatrix>{xf})[0].toScalar();
//    result.Q = result.Ef - result.Ei;
//
//
//    vector<double> pi(L);
//    result.p = 0;
//    for (int i = 0; i < L; i++) {
//        pi[i] = 1 - norm(dot(ff[i], f0[i]));
//        result.p += pi[i];
//    }
//    result.p /= L;
//
//    ptime stop_time = microsec_clock::local_time();
//    time_period period(start_time, stop_time);
//    result.runtime = to_simple_string(period.length());
//}

SX DynamicsProblem::energy(SX& fin, SX& J, SX& U0, SX& dU, double mu) {

    SX E = 0;
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            E += energy(i, n, fin, J, U0, dU, mu);
        }
    }
    return E;
}

SX DynamicsProblem::energy(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, double mu) {

    complex<SX> expth = complex<SX>(1, 0);
    complex<SX> expmth = ~expth;
    complex<SX> exp2th = expth*expth;
    complex<SX> expm2th = ~exp2th;

    vector<vector<complex < SX>>> f(L, vector<complex < SX >> (dim, complex<SX>(0, 0)));
    vector<SX> norm2(L, 0);
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < dim; k++) {
            int l = 2 * (j * dim + k);
            f[j][k] = complex<SX>(fin[l], fin[l + 1]);
        }
        for (int m = 0; m <= nmax; m++) {
            norm2[j] += f[j][m].real() * f[j][m].real() + f[j][m].imag() * f[j][m].imag();
        }
    }

    complex<SX> E = complex<SX>(0, 0);

    complex<SX> Ei, Ej1, Ej2, Ej1j2, Ej1k1, Ej2k2;

    int k1 = mod(i - 2);
    int j1 = mod(i - 1);
    int j2 = mod(i + 1);
    int k2 = mod(i + 2);

    Ei = complex<SX>(0, 0);
    Ej1 = complex<SX>(0, 0);
    Ej2 = complex<SX>(0, 0);
    Ej1j2 = complex<SX>(0, 0);
    Ej1k1 = complex<SX>(0, 0);
    Ej2k2 = complex<SX>(0, 0);

    Ei += (0.5 * (U0 + dU[i]) * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

    if (n < nmax) {
        Ej1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                * f[i][n] * f[j1][n + 1];
        Ej2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
                * f[j2][n + 1];

        if (n > 0) {
            Ej1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j1][n - 1] * f[i][n - 1] * f[j1][n + 1];
            Ej2 += 0.5 * J[i] * J[i] * expm2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j2][n - 1] * f[i][n - 1] * f[j2][n + 1];
        }
        if (n < nmax - 1) {
            Ej1 -= 0.5 * J[j1] * J[j1] * exp2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j1][n] * f[i][n] * f[j1][n + 2];
            Ej2 -= 0.5 * J[i] * J[i] * expm2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j2][n] * f[i][n] * f[j2][n + 2];
        }

        if (n > 1) {
            Ej1 += -J[j1] * J[j1] * exp2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j1, n, n - 1, i, j1, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j1][n - 2] * f[i][n - 1] * f[j1][n];
            Ej2 += -J[i] * J[i] * expm2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j2, n, n - 1, i, j2, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j2][n - 2] * f[i][n - 1] * f[j2][n];
        }
        if (n < nmax - 2) {
            Ej1 -= -J[j1] * J[j1] * exp2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j1, n, n + 3, i, j1, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j1][n + 1] * f[i][n] * f[j1][n + 3];
            Ej2 -= -J[i] * J[i] * expm2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j2, n, n + 3, i, j2, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j2][n + 1] * f[i][n] * f[j2][n + 3];
        }

        for (int m = 1; m <= nmax; m++) {
            if (n != m - 1) {
                Ej1 += 0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1] -
                        ~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m]);
                Ej2 += 0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1] -
                        ~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m]);

                Ej1 += 1.0 * J[j1] * expth * g(n, m) * (eps(dU, i, j1, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n] * f[j1][m];
                Ej2 += 1.0 * J[i] * expmth * g(n, m) * (eps(dU, i, j2, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n] * f[j2][m];

                if (n != m - 3 && m > 1 && n < nmax - 1) {
                    Ej1 += -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 += -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
                if (n != m + 1 && n > 0 && m < nmax) {
                    Ej1 -= -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 -= -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }

                if (n > 0) {
                    Ej1j2 += -J[j1] * J[i] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j1, n, m, i, j2, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][n - 1]
                            * f[i][n - 1] * f[j1][m] * f[j2][n];
                    Ej1j2 += -J[i] * J[j1] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j2, n, m, i, j1, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][n - 1]
                            * f[i][n - 1] * f[j2][m] * f[j1][n];
                }
                if (n < nmax - 1) {
                    Ej1j2 -= -J[j1] * J[i] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j1, n, m, i, j2, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][n + 1]
                            * f[i][n] * f[j1][m] * f[j2][n + 2];
                    Ej1j2 -= -J[i] * J[j1] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j2, n, m, i, j1, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][n + 1]
                            * f[i][n] * f[j2][m] * f[j1][n + 2];
                }

                Ej1 += -0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m] -
                        ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1]);
                Ej2 += -0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m] -
                        ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1]);

                for (int q = 1; q <= nmax; q++) {
                    if (n < nmax - 1 && n != q - 2) {
                        Ej1j2 += -0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n] * f[j1][m] * f[j2][q];
                        Ej1j2 += -0.5 * J[i] * J[j1] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n] * f[j2][m] * f[j1][q];
                    }
                    if (n > 0 && n != q) {
                        Ej1j2 -= -0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n - 1] * f[j1][m] * f[j2][q];
                        Ej1j2 -= -0.5 * J[i] * J[j1] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n - 1] * f[j2][m] * f[j1][q];
                    }

                    if (m != q) {
                        Ej1k1 += -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m] * f[k1][q];
                        Ej2k2 += -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m] * f[k2][q];
                        Ej1k1 -= -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m - 1] * f[k1][q];
                        Ej2k2 -= -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m - 1] * f[k2][q];
                    }

                }

                for (int p = 0; p < nmax; p++) {

                    if (p != n - 1 && 2 * n - m == p && n > 0) {
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m == p - 2 && n < nmax - 1) {
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != n - 1 && 2 * n - m != p && n > 0) {
                        Ej1j2 += -0.25 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += -0.25 * J[i] * J[j1] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m != p - 2 && n < nmax - 1) {
                        Ej1j2 -= -0.25 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= -0.25 * J[i] * J[j1] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != m - 1 && n != p) {
                        Ej1k1 += -0.25 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j1, n, m, j1, k1, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][p] * f[i][n] * f[j1][m - 1] * f[k1][p + 1] -
                                ~f[i][n + 1] * ~f[j1][m] * ~f[k1][p] * f[i][n] * f[j1][m] * f[k1][p + 1]);
                        Ej2k2 += -0.25 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j2, n, m, j2, k2, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][p] * f[i][n] * f[j2][m - 1] * f[k2][p + 1] -
                                ~f[i][n + 1] * ~f[j2][m] * ~f[k2][p] * f[i][n] * f[j2][m] * f[k2][p + 1]);
                    }
                }

                Ej1k1 += 0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][n]
                        * f[i][n] * f[j1][m - 1] * f[k1][n + 1] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][n]
                        * f[i][n] * f[j1][m] * f[k1][n + 1]);
                Ej2k2 += 0.5 * J[j2] * J[i] * expm2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][n]
                        * f[i][n] * f[j2][m - 1] * f[k2][n + 1] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][n]
                        * f[i][n] * f[j2][m] * f[k2][n + 1]);

                Ej1k1 += -J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j1, n, m, j1, k1, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][m - 1] * f[i][n] * f[j1][m - 1] * f[k1][m] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][m - 1] * f[i][n] * f[j1][m] * f[k1][m]);
                Ej2k2 += -J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j2, n, m, j2, k2, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][m - 1] * f[i][n] * f[j2][m - 1] * f[k2][m] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][m - 1] * f[i][n] * f[j2][m] * f[k2][m]);

                if (m != n - 1 && n != m && m < nmax && n > 0) {
                    Ej1 += -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m, i, j1, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 += -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m, i, j2, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }
                if (n != m - 3 && n != m - 2 && n < nmax - 1 && m > 1) {
                    Ej1 -= -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m, i, j1, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 -= -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m, i, j2, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
            }
        }
    }

    Ei /= norm2[i];
    Ej1 /= norm2[i] * norm2[j1];
    Ej2 /= norm2[i] * norm2[j2];
    Ej1j2 /= norm2[i] * norm2[j1] * norm2[j2];
    Ej1k1 /= norm2[i] * norm2[j1] * norm2[k1];
    Ej2k2 /= norm2[i] * norm2[j2] * norm2[k2];

    E += Ei;
    E += Ej1;
    E += Ej2;
    E += Ej1j2;
    E += Ej1k1;
    E += Ej2k2;

    return E.real();
}

SX DynamicsProblem::canonical(SX& fin, SX& J, SX& U0, SX& dU, double mu) {

    SX S = 0;
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            S += canonical(i, n, fin, J, U0, dU, mu);
        }
    }
    return S;
}

SX DynamicsProblem::canonical(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, double mu) {

    vector<vector<complex < SX>>> f(L, vector<complex < SX >> (dim, complex<SX>(0, 0)));
    vector<SX> norm2(L, 0);
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < dim; k++) {
            int l = 2 * (j * dim + k);
            f[j][k] = complex<SX>(fin[l], fin[l + 1]);
        }
        for (int m = 0; m <= nmax; m++) {
            norm2[j] += f[j][m].real() * f[j][m].real() + f[j][m].imag() * f[j][m].imag();
        }
    }

    complex<SX> S = complex<SX>(0, 0);

    complex<SX> Sj1, Sj2;

    int j1 = mod(i - 1);
    int j2 = mod(i + 1);

    Sj1 = complex<SX>(0, 0);
    Sj2 = complex<SX>(0, 0);

    if (n < nmax) {
        for (int m = 1; m <= nmax; m++) {
            if (n != m - 1) {
                Sj1 += -J[j1] * g(n, m) * (1 / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n] * f[j1][m];
                Sj2 += -J[i] * g(n, m) * (1 / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n] * f[j2][m];

            }
        }
    }

    Sj1 /= norm2[i] * norm2[j1];
    Sj2 /= norm2[i] * norm2[j2];

    S += Sj1;
    S += Sj2;

    return S.imag();
}

template<class T> SX energy(SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    SX E = 0;
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            E += energy(i, n, fin, J, U0, dU, mu);
        }
    }
    return E;
}

template<class T> SX energy(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    complex<SX> expth = complex<SX>(1, 0);
    complex<SX> expmth = ~expth;
    complex<SX> exp2th = expth*expth;
    complex<SX> expm2th = ~exp2th;

    vector<vector<complex < SX>>> f(L, vector<complex < SX >> (dim, complex<SX>(0, 0)));
    vector<SX> norm2(L, 0);
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < dim; k++) {
            int l = 2 * (j * dim + k);
            f[j][k] = complex<SX>(fin[l], fin[l + 1]);
        }
        for (int m = 0; m <= nmax; m++) {
            norm2[j] += f[j][m].real() * f[j][m].real() + f[j][m].imag() * f[j][m].imag();
        }
    }

    complex<SX> E = complex<SX>(0, 0);

    complex<SX> Ei, Ej1, Ej2, Ej1j2, Ej1k1, Ej2k2;

    int k1 = mod(i - 2);
    int j1 = mod(i - 1);
    int j2 = mod(i + 1);
    int k2 = mod(i + 2);

    Ei = complex<SX>(0, 0);
    Ej1 = complex<SX>(0, 0);
    Ej2 = complex<SX>(0, 0);
    Ej1j2 = complex<SX>(0, 0);
    Ej1k1 = complex<SX>(0, 0);
    Ej2k2 = complex<SX>(0, 0);

    Ei += (0.5 * (U0 + dU[i]) * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

    if (n < nmax) {
        Ej1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                * f[i][n] * f[j1][n + 1];
        Ej2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
                * f[j2][n + 1];

        if (n > 0) {
            Ej1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j1][n - 1] * f[i][n - 1] * f[j1][n + 1];
            Ej2 += 0.5 * J[i] * J[i] * expm2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j2][n - 1] * f[i][n - 1] * f[j2][n + 1];
        }
        if (n < nmax - 1) {
            Ej1 -= 0.5 * J[j1] * J[j1] * exp2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j1][n] * f[i][n] * f[j1][n + 2];
            Ej2 -= 0.5 * J[i] * J[i] * expm2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j2][n] * f[i][n] * f[j2][n + 2];
        }

        if (n > 1) {
            Ej1 += -J[j1] * J[j1] * exp2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j1, n, n - 1, i, j1, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j1][n - 2] * f[i][n - 1] * f[j1][n];
            Ej2 += -J[i] * J[i] * expm2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j2, n, n - 1, i, j2, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j2][n - 2] * f[i][n - 1] * f[j2][n];
        }
        if (n < nmax - 2) {
            Ej1 -= -J[j1] * J[j1] * exp2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j1, n, n + 3, i, j1, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j1][n + 1] * f[i][n] * f[j1][n + 3];
            Ej2 -= -J[i] * J[i] * expm2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j2, n, n + 3, i, j2, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j2][n + 1] * f[i][n] * f[j2][n + 3];
        }

        for (int m = 1; m <= nmax; m++) {
            if (n != m - 1) {
                Ej1 += 0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1] -
                        ~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m]);
                Ej2 += 0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1] -
                        ~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m]);

                Ej1 += 1.0 * J[j1] * expth * g(n, m) * (eps(dU, i, j1, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n] * f[j1][m];
                Ej2 += 1.0 * J[i] * expmth * g(n, m) * (eps(dU, i, j2, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n] * f[j2][m];

                if (n != m - 3 && m > 1 && n < nmax - 1) {
                    Ej1 += -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 += -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
                if (n != m + 1 && n > 0 && m < nmax) {
                    Ej1 -= -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 -= -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }

                if (n > 0) {
                    Ej1j2 += -J[j1] * J[i] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j1, n, m, i, j2, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][n - 1]
                            * f[i][n - 1] * f[j1][m] * f[j2][n];
                    Ej1j2 += -J[i] * J[j1] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j2, n, m, i, j1, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][n - 1]
                            * f[i][n - 1] * f[j2][m] * f[j1][n];
                }
                if (n < nmax - 1) {
                    Ej1j2 -= -J[j1] * J[i] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j1, n, m, i, j2, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][n + 1]
                            * f[i][n] * f[j1][m] * f[j2][n + 2];
                    Ej1j2 -= -J[i] * J[j1] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j2, n, m, i, j1, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][n + 1]
                            * f[i][n] * f[j2][m] * f[j1][n + 2];
                }

                Ej1 += -0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m] -
                        ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1]);
                Ej2 += -0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m] -
                        ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1]);

                for (int q = 1; q <= nmax; q++) {
                    if (n < nmax - 1 && n != q - 2) {
                        Ej1j2 += -0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n] * f[j1][m] * f[j2][q];
                        Ej1j2 += -0.5 * J[i] * J[j1] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n] * f[j2][m] * f[j1][q];
                    }
                    if (n > 0 && n != q) {
                        Ej1j2 -= -0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n - 1] * f[j1][m] * f[j2][q];
                        Ej1j2 -= -0.5 * J[i] * J[j1] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n - 1] * f[j2][m] * f[j1][q];
                    }

                    if (m != q) {
                        Ej1k1 += -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m] * f[k1][q];
                        Ej2k2 += -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m] * f[k2][q];
                        Ej1k1 -= -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m - 1] * f[k1][q];
                        Ej2k2 -= -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m - 1] * f[k2][q];
                    }

                }

                for (int p = 0; p < nmax; p++) {

                    if (p != n - 1 && 2 * n - m == p && n > 0) {
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m == p - 2 && n < nmax - 1) {
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != n - 1 && 2 * n - m != p && n > 0) {
                        Ej1j2 += -0.25 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += -0.25 * J[i] * J[j1] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m != p - 2 && n < nmax - 1) {
                        Ej1j2 -= -0.25 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= -0.25 * J[i] * J[j1] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != m - 1 && n != p) {
                        Ej1k1 += -0.25 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j1, n, m, j1, k1, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][p] * f[i][n] * f[j1][m - 1] * f[k1][p + 1] -
                                ~f[i][n + 1] * ~f[j1][m] * ~f[k1][p] * f[i][n] * f[j1][m] * f[k1][p + 1]);
                        Ej2k2 += -0.25 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j2, n, m, j2, k2, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][p] * f[i][n] * f[j2][m - 1] * f[k2][p + 1] -
                                ~f[i][n + 1] * ~f[j2][m] * ~f[k2][p] * f[i][n] * f[j2][m] * f[k2][p + 1]);
                    }
                }

                Ej1k1 += 0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][n]
                        * f[i][n] * f[j1][m - 1] * f[k1][n + 1] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][n]
                        * f[i][n] * f[j1][m] * f[k1][n + 1]);
                Ej2k2 += 0.5 * J[j2] * J[i] * expm2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][n]
                        * f[i][n] * f[j2][m - 1] * f[k2][n + 1] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][n]
                        * f[i][n] * f[j2][m] * f[k2][n + 1]);

                Ej1k1 += -J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j1, n, m, j1, k1, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][m - 1] * f[i][n] * f[j1][m - 1] * f[k1][m] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][m - 1] * f[i][n] * f[j1][m] * f[k1][m]);
                Ej2k2 += -J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j2, n, m, j2, k2, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][m - 1] * f[i][n] * f[j2][m - 1] * f[k2][m] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][m - 1] * f[i][n] * f[j2][m] * f[k2][m]);

                if (m != n - 1 && n != m && m < nmax && n > 0) {
                    Ej1 += -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m, i, j1, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 += -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m, i, j2, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }
                if (n != m - 3 && n != m - 2 && n < nmax - 1 && m > 1) {
                    Ej1 -= -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m, i, j1, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 -= -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m, i, j2, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
            }
        }
    }

    Ei /= norm2[i];
    Ej1 /= norm2[i] * norm2[j1];
    Ej2 /= norm2[i] * norm2[j2];
    Ej1j2 /= norm2[i] * norm2[j1] * norm2[j2];
    Ej1k1 /= norm2[i] * norm2[j1] * norm2[k1];
    Ej2k2 /= norm2[i] * norm2[j2] * norm2[k2];

    E += Ei;
    E += Ej1;
    E += Ej2;
    E += Ej1j2;
    E += Ej1k1;
    E += Ej2k2;

    return E.real();
}

template<class T> SX canonical(SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    SX S = 0;
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            S += canonical(i, n, fin, J, U0, dU, mu);
        }
    }
    return S;
}

template<class T> SX canonical(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    vector<vector<complex < SX>>> f(L, vector<complex < SX >> (dim, complex<SX>(0, 0)));
    vector<SX> norm2(L, 0);
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < dim; k++) {
            int l = 2 * (j * dim + k);
            f[j][k] = complex<SX>(fin[l], fin[l + 1]);
        }
        for (int m = 0; m <= nmax; m++) {
            norm2[j] += f[j][m].real() * f[j][m].real() + f[j][m].imag() * f[j][m].imag();
        }
    }

    complex<SX> S = complex<SX>(0, 0);

    complex<SX> Sj1, Sj2;

    int j1 = mod(i - 1);
    int j2 = mod(i + 1);

    Sj1 = complex<SX>(0, 0);
    Sj2 = complex<SX>(0, 0);

    if (n < nmax) {
        for (int m = 1; m <= nmax; m++) {
            if (n != m - 1) {
                Sj1 += -J[j1] * g(n, m) * (1 / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n] * f[j1][m];
                Sj2 += -J[i] * g(n, m) * (1 / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n] * f[j2][m];

            }
        }
    }

    Sj1 /= norm2[i] * norm2[j1];
    Sj2 /= norm2[i] * norm2[j2];

    S += Sj1;
    S += Sj2;

    return S.imag();
}

template<class T> SX energy2(SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    SX E = 0;
    for (int i = 0; i < L; i++) {
        for (int n = 0; n <= nmax; n++) {
            E += energy2(i, n, fin, J, U0, dU, mu);
        }
    }
    return E;
}

template<class T> SX energy2(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, T mu) {

    complex<SX> expth = complex<SX>(1, 0);
//    complex<SX> expth = complex<SX>(sin(0.1), cos(0.1));
    complex<SX> expmth = ~expth;
    complex<SX> exp2th = expth*expth;
    complex<SX> expm2th = ~exp2th;

    vector<vector<complex < SX>>> f(L, vector<complex < SX >> (dim, complex<SX>(0, 0)));
    vector<SX> norm2(L, 0);
    for (int j = 0; j < L; j++) {
        for (int k = 0; k < dim; k++) {
            int l = 2 * (j * dim + k);
            f[j][k] = complex<SX>(fin[l], fin[l + 1]);
        }
        for (int m = 0; m <= nmax; m++) {
            norm2[j] += f[j][m].real() * f[j][m].real() + f[j][m].imag() * f[j][m].imag();
        }
    }

    complex<SX> E = complex<SX>(0, 0);

    complex<SX> Ei, Ej1, Ej2, Ej1j2, Ej1k1, Ej2k2;

    int k1 = mod(i - 2);
    int j1 = mod(i - 1);
    int j2 = mod(i + 1);
    int k2 = mod(i + 2);

    Ei = complex<SX>(0, 0);
    Ej1 = complex<SX>(0, 0);
    Ej2 = complex<SX>(0, 0);
    Ej1j2 = complex<SX>(0, 0);
    Ej1k1 = complex<SX>(0, 0);
    Ej2k2 = complex<SX>(0, 0);

    Ei += (0.5 * (U0 + dU[i]) * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

    if (n < nmax) {
        Ej1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                * f[i][n] * f[j1][n + 1];
        Ej2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
                * f[j2][n + 1];

        if (n > 0) {
            Ej1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j1][n - 1] * f[i][n - 1] * f[j1][n + 1];
            Ej2 += 0.5 * J[i] * J[i] * expm2th * g(n, n) * g(n - 1, n + 1) * (1 / eps(U0, n, n))
                    * ~f[i][n + 1] * ~f[j2][n - 1] * f[i][n - 1] * f[j2][n + 1];
        }
        if (n < nmax - 1) {
            Ej1 -= 0.5 * J[j1] * J[j1] * exp2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j1][n] * f[i][n] * f[j1][n + 2];
            Ej2 -= 0.5 * J[i] * J[i] * expm2th * g(n, n + 2) * g(n + 1, n + 1) * (1 / eps(U0, n, n + 2))
                    * ~f[i][n + 2] * ~f[j2][n] * f[i][n] * f[j2][n + 2];
        }

        if (n > 1) {
            Ej1 += -J[j1] * J[j1] * exp2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j1, n, n - 1, i, j1, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j1][n - 2] * f[i][n - 1] * f[j1][n];
            Ej2 += -J[i] * J[i] * expm2th * g(n, n - 1) * g(n - 1, n)
                    * (eps(dU, i, j2, n, n - 1, i, j2, n - 1, n) / (eps(U0, n, n - 1)*(eps(U0, n, n - 1) + eps(U0, n - 1, n))))
                    * ~f[i][n + 1] * ~f[j2][n - 2] * f[i][n - 1] * f[j2][n];
        }
        if (n < nmax - 2) {
            Ej1 -= -J[j1] * J[j1] * exp2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j1, n, n + 3, i, j1, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j1][n + 1] * f[i][n] * f[j1][n + 3];
            Ej2 -= -J[i] * J[i] * expm2th * g(n, n + 3) * g(n + 1, n + 2)
                    * (eps(dU, i, j2, n, n + 3, i, j2, n + 1, n + 2) / (eps(U0, n, n + 3)*(eps(U0, n, n + 3) + eps(U0, n + 1, n + 2))))
                    * ~f[i][n + 2] * ~f[j2][n + 1] * f[i][n] * f[j2][n + 3];
        }

        for (int m = 1; m <= nmax; m++) {
            if (n != m - 1) {
                Ej1 += 0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1] -
                        ~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m]);
                Ej2 += 0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1) * (1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1] -
                        ~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m]);

                Ej1 += 1.0 * J[j1] * expth * g(n, m) * (eps(dU, i, j1, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n] * f[j1][m];
                Ej2 += 1.0 * J[i] * expmth * g(n, m) * (eps(dU, i, j2, n, m) / eps(U0, n, m))
                        * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n] * f[j2][m];

                if (n != m - 3 && m > 1 && n < nmax - 1) {
                    Ej1 += -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 += -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, m - 1)))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
                if (n != m + 1 && n > 0 && m < nmax) {
                    Ej1 -= -0.5 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 -= -0.5 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, m + 1)))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }

                if (n > 0) {
                    Ej1j2 += -J[j1] * J[i] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j1, n, m, i, j2, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][n - 1]
                            * f[i][n - 1] * f[j1][m] * f[j2][n];
                    Ej1j2 += -J[i] * J[j1] * g(n, m) * g(n - 1, n)
                            * (eps(dU, i, j2, n, m, i, j1, n - 1, n) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n - 1, n))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][n - 1]
                            * f[i][n - 1] * f[j2][m] * f[j1][n];
                }
                if (n < nmax - 1) {
                    Ej1j2 -= -J[j1] * J[i] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j1, n, m, i, j2, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][n + 1]
                            * f[i][n] * f[j1][m] * f[j2][n + 2];
                    Ej1j2 -= -J[i] * J[j1] * g(n, m) * g(n + 1, n + 2)
                            * (eps(dU, i, j2, n, m, i, j1, n + 1, n + 2) / (eps(U0, n, m) * (eps(U0, n, m) + eps(U0, n + 1, n + 2))))
                            * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][n + 1]
                            * f[i][n] * f[j2][m] * f[j1][n + 2];
                }

                Ej1 += -0.5 * J[j1] * J[j1] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m] -
                        ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1]);
                Ej2 += -0.5 * J[i] * J[i] * g(n, m) * g(m - 1, n + 1)
                        * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, n + 1)))
                        * (~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m] -
                        ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1]);

                for (int q = 1; q <= nmax; q++) {
                    if (n < nmax - 1 && n != q - 2) {
                        Ej1j2 += -0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n] * f[j1][m] * f[j2][q];
                        Ej1j2 += -0.5 * J[i] * J[j1] * g(n, m) * g(n + 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n + 1, q)))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n] * f[j2][m] * f[j1][q];
                    }
                    if (n > 0 && n != q) {
                        Ej1j2 -= -0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][q - 1]
                                * f[i][n - 1] * f[j1][m] * f[j2][q];
                        Ej1j2 -= -0.5 * J[i] * J[j1] * g(n, m) * g(n - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, n - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][q - 1]
                                * f[i][n - 1] * f[j2][m] * f[j1][q];
                    }

                    if (m != q) {
                        Ej1k1 += -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m] * f[k1][q];
                        Ej2k2 += -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m] * f[k2][q];
                        Ej1k1 -= -0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j1, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][q - 1]
                                * f[i][n] * f[j1][m - 1] * f[k1][q];
                        Ej2k2 -= -0.5 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, q)
                                * (eps(dU, i, j2, n, m) / (eps(U0, n, m) * eps(U0, m - 1, q)))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][q - 1]
                                * f[i][n] * f[j2][m - 1] * f[k2][q];
                    }

                }

                for (int p = 0; p < nmax; p++) {

                    if (p != n - 1 && 2 * n - m == p && n > 0) {
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += 0.5 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m == p - 2 && n < nmax - 1) {
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= 0.5 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1) * (1 / eps(U0, n, m))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != n - 1 && 2 * n - m != p && n > 0) {
                        Ej1j2 += -0.25 * J[j1] * J[i] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n - 1] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 += -0.25 * J[i] * J[j1] * g(n, m) * g(n - 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, p + 1))))
                                * ~f[i][n + 1] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n - 1] * f[j2][m] * f[j1][p + 1];
                    }
                    if (p != n + 1 && 2 * n - m != p - 2 && n < nmax - 1) {
                        Ej1j2 -= -0.25 * J[j1] * J[i] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j1, n, m, i, j2, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j1][m - 1] * ~f[j2][p]
                                * f[i][n] * f[j1][m] * f[j2][p + 1];
                        Ej1j2 -= -0.25 * J[i] * J[j1] * g(n, m) * g(n + 1, p + 1)
                                * (eps(dU, i, j2, n, m, i, j1, p, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, p + 1))))
                                * ~f[i][n + 2] * ~f[j2][m - 1] * ~f[j1][p]
                                * f[i][n] * f[j2][m] * f[j1][p + 1];
                    }

                    if (p != m - 1 && n != p) {
                        Ej1k1 += -0.25 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j1, n, m, j1, k1, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][p] * f[i][n] * f[j1][m - 1] * f[k1][p + 1] -
                                ~f[i][n + 1] * ~f[j1][m] * ~f[k1][p] * f[i][n] * f[j1][m] * f[k1][p + 1]);
                        Ej2k2 += -0.25 * J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, p + 1)
                                * (eps(dU, i, j2, n, m, j2, k2, p, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, p + 1))))
                                * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][p] * f[i][n] * f[j2][m - 1] * f[k2][p + 1] -
                                ~f[i][n + 1] * ~f[j2][m] * ~f[k2][p] * f[i][n] * f[j2][m] * f[k2][p + 1]);
                    }
                }

                Ej1k1 += 0.5 * J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][n]
                        * f[i][n] * f[j1][m - 1] * f[k1][n + 1] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][n]
                        * f[i][n] * f[j1][m] * f[k1][n + 1]);
                Ej2k2 += 0.5 * J[j2] * J[i] * expm2th * g(n, m) * g(m - 1, n + 1)*(1 / eps(U0, n, m))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][n]
                        * f[i][n] * f[j2][m - 1] * f[k2][n + 1] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][n]
                        * f[i][n] * f[j2][m] * f[k2][n + 1]);

                Ej1k1 += -J[j1] * J[k1] * exp2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j1, n, m, j1, k1, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j1][m - 1] * ~f[k1][m - 1] * f[i][n] * f[j1][m - 1] * f[k1][m] -
                        ~f[i][n + 1] * ~f[j1][m] * ~f[k1][m - 1] * f[i][n] * f[j1][m] * f[k1][m]);
                Ej2k2 += -J[i] * J[j2] * expm2th * g(n, m) * g(m - 1, m)
                        * (eps(dU, i, j2, n, m, j2, k2, m - 1, m) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, m - 1, m))))
                        * (~f[i][n + 1] * ~f[j2][m - 1] * ~f[k2][m - 1] * f[i][n] * f[j2][m - 1] * f[k2][m] -
                        ~f[i][n + 1] * ~f[j2][m] * ~f[k2][m - 1] * f[i][n] * f[j2][m] * f[k2][m]);

                if (m != n - 1 && n != m && m < nmax && n > 0) {
                    Ej1 += -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j1, n, m, i, j1, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j1][m - 1] * f[i][n - 1] * f[j1][m + 1];
                    Ej2 += -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n - 1, m + 1)
                            * (eps(dU, i, j2, n, m, i, j2, m, n) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n - 1, m + 1))))
                            * ~f[i][n + 1] * ~f[j2][m - 1] * f[i][n - 1] * f[j2][m + 1];
                }
                if (n != m - 3 && n != m - 2 && n < nmax - 1 && m > 1) {
                    Ej1 -= -0.25 * J[j1] * J[j1] * exp2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j1, n, m, i, j1, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j1][m - 2] * f[i][n] * f[j1][m];
                    Ej2 -= -0.25 * J[i] * J[i] * expm2th * g(n, m) * g(n + 1, m - 1)
                            * (eps(dU, i, j2, n, m, i, j2, m - 2, n + 2) / (eps(U0, n, m)*(eps(U0, n, m) + eps(U0, n + 1, m - 1))))
                            * ~f[i][n + 2] * ~f[j2][m - 2] * f[i][n] * f[j2][m];
                }
            }
        }
    }

    Ei /= norm2[i];
    Ej1 /= norm2[i] * norm2[j1];
    Ej2 /= norm2[i] * norm2[j2];
    Ej1j2 /= norm2[i] * norm2[j1] * norm2[j2];
    Ej1k1 /= norm2[i] * norm2[j1] * norm2[k1];
    Ej2k2 /= norm2[i] * norm2[j2] * norm2[k2];

    E += Ei;
    E += Ej1;
    E += Ej2;
    E += Ej1j2;
    E += Ej1k1;
    E += Ej2k2;

    return E.real();
}

