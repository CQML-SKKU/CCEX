/**
 * @file simulator_cce_drive.cpp
 * @brief cCCE with continuous bath driving (pulse.bath_drive): rotating frames, rotating-wave approximation.
 *
 * Port of pyCCEX ccex/simulator/drive.py (feature/spin-squeezing b3ce43e), ensemble approach (no s_range).
 * docs/methods/bath-drive.md has the model. In short: every bath spin k sits on a line with lab frequency
 * w0 + d_k; a tone drives a set of lines at w0 + line + offset. A driven spin is put in the frame of its tone, the
 * others in the frame of their own line, and only what is static there is kept:
 *     H_k  = (d_k - f_k) Jz_k + W (cos p Jx_k + sin p Jy_k)
 *     H_kj = T_zz [Jz_k Jz_j - (Jx_k Jx_j + Jy_k Jy_j)/2]   same frame,   T_zz Jz_k Jz_j   different frames
 *     H_qb = m_s sum_k A_zz,k Jz_k
 * With qubit_transverse (default) the qubit coupling A_zx, A_zy is kept exactly: in branch m_s spin k precesses about
 * w_k z + m_s A_k, each branch gets its own tilted axes, and the frame changes at the pi pulses and at the end are
 * applied. shifts adds the leading lab corrections (Bloch-Siegert, AC Stark, Van Vleck, drive-assisted cross-line
 * flip-flop). Units: rad*kHz, ms, Angstrom, G.
 */
#include "../../include/calculation/hamiltonian.h"
#include "../../include/calculation/simulator.h"
#include "../../include/calculation/bath_drive.h"
#include "../../include/input/reader.h"
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr double kFrameEps = 1e-6;   // rad*kHz: two frames are the same only when built from the same line + offset

struct Component { int tone; double line, carrier, W, phase; };   // one (tone, listed line)
struct Frame { double f, dz, wx, wy; };                             // frame, static Jz coefficient, field x, y
struct Axes {                                                       // one qubit branch, per bath spin
    std::vector<Eigen::Vector3d> n, e1, e2;
    std::vector<double> delta;                                      // W_k - w_k along n_k
};

std::set<std::string>& warned(){ static std::set<std::string> s; return s; }
void warnOnce(const std::string& key, const std::string& msg){
    if (warned().insert(key).second) fprintf(stderr, "Warning (rank %d): bath drive: %s\n", rank, msg.c_str());
}

double fieldZ(Config* cnf){ return (double)Config_getBfield(cnf)[2]; }

bool sameComponent(const Component& a, const Component& b){ return a.tone == b.tone && a.line == b.line; }

std::vector<Component> components(const BathDrive& dr){
    std::vector<Component> out;
    for (size_t it = 0; it < dr.tones.size(); ++it){
        const BathDriveTone& t = dr.tones[it];
        if (t.rabi == 0.0) continue;                    // no field: must not move any spin into a frame
        std::vector<double> seen;
        for (double l : t.lines){
            if (std::find(seen.begin(), seen.end(), l) != seen.end()) continue;   // a repeated line counts once
            seen.push_back(l);
            out.push_back({(int)it, l, l + t.offset, t.rabi, t.phase});
        }
    }
    return out;
}

struct SpinData { double gyro, d, spin; std::string name; };     // what the frames need of a spin

std::vector<SpinData> spinData(BathArray* ba){
    std::vector<SpinData> out(BathArray_getNspin(ba));
    for (size_t k = 0; k < out.size(); ++k){
        BathSpin* bs = ba->bath[k];
        out[k] = {BathSpin_getGyro(bs), BathSpin_getDetuning(bs), (double)BathSpin_getSpin(bs), BathSpin_getName(bs)};
    }
    return out;
}

// Per spin: (frame f, static Jz coefficient, Wx, Wy). check: the refusals and warnings on the whole bath
std::vector<Frame> frameAssignment(const std::vector<SpinData>& sp, Config* cnf, const BathDrive& dr,
                                   const std::vector<Component>& comps, bool check){
    const double bz = fieldZ(cnf);
    const double tol = dr.line_tolerance;
    const int n = (int)sp.size();
    std::vector<Frame> out(n);
    for (int k = 0; k < n; ++k){
        const double gyro = sp[k].gyro, d = sp[k].d;
        if (check){
            if (std::fabs(sp[k].spin - 0.5) > 1e-9 || gyro >= 0.0){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive drives electron-like spin-1/2 bath spins only (gyro < 0); "
                                 "spin %d (%s) has spin %g, gyro %g\n", k, sp[k].name.c_str(), sp[k].spin, gyro);
                ccexRefuse();
            }
            if (!std::isfinite(d) || !std::isfinite(gyro)){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive: bath spin %d has a non-finite detuning or gyro\n", k);
                ccexRefuse();
            }
        }
        const double w0 = -gyro * bz;                   // electron Zeeman frequency, > 0
        const double wk = w0 + d;
        if (check) for (const Component& c : comps)
            if (w0 + c.carrier <= 3.0 * c.W || wk <= 3.0 * c.W){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive: a lab carrier or line frequency is within 3 W of zero "
                                 "(counter-rotating resonance); the rotating-frame model does not apply\n");
                ccexRefuse();
            }
        const Component* mine = nullptr;
        std::set<int> ownTones;
        for (const Component& c : comps)
            if (std::fabs(d - c.line) < tol){
                ownTones.insert(c.tone);
                if (!mine || std::fabs(d - c.line) < std::fabs(d - mine->line)) mine = &c;
            }
        if (check && ownTones.size() > 1){
            ccexRefusalPrint(stderr, "Error: pulse.bath_drive: bath spin %d (detuning %.6f MHz) is driven by more than one tone\n",
                             k, d / MHZ_TO_RADKHZ(1.0));
            ccexRefuse();
        }
        const double Wown = mine ? mine->W : 0.0;
        if (check) for (const Component& c : comps)
            if ((!mine || !sameComponent(c, *mine)) && c.W > 0.0){
                const double gap = d - c.carrier;
                if (std::fabs(gap) <= 3.0 * std::max(c.W, Wown)){
                    char key[64], msg[256];
                    snprintf(key, sizeof key, "near %.3f", gap);
                    snprintf(msg, sizeof msg, "a drive carrier is within 3 W of a line it does not target (gap %.3f MHz); "
                             "the per-line RWA drops a near-resonant process", gap / MHZ_TO_RADKHZ(1.0));
                    warnOnce(key, msg);
                }
            }
        if (!mine){ out[k] = {d, 0.0, 0.0, 0.0}; continue; }
        double dz = d - mine->carrier;
        if (dr.shifts) for (const Component& c : comps){
            const double wc = w0 + c.carrier;
            dz += c.W * c.W / (2.0 * (wk + wc));                              // Bloch-Siegert
            if (!sameComponent(c, *mine) && std::fabs(wk - wc) > 3.0 * c.W)
                dz += c.W * c.W / (2.0 * (wk - wc));                          // AC Stark of the other components
        }
        out[k] = {mine->carrier, dz, mine->W * std::cos(mine->phase), mine->W * std::sin(mine->phase)};
    }
    return out;
}

std::vector<Frame> frameAssignment(BathArray* ba, Config* cnf, const BathDrive& dr, const std::vector<Component>& comps,
                                   bool check){
    return frameAssignment(spinData(ba), cnf, dr, comps, check);
}

// ---- operators on the 2^n space of a cluster: spin 0 is the first Kronecker factor (the highest bit) ----------------

inline int bitOf(int b, int k, int n){ return (b >> (n - 1 - k)) & 1; }

// H += op (2x2) on spin k
void addOne(MatrixXcd& H, const Eigen::Matrix2cd& op, int k, int n){
    const int d = 1 << n, mask = 1 << (n - 1 - k);
    for (int b = 0; b < d; ++b){
        const int s = bitOf(b, k, n);
        for (int r = 0; r < 2; ++r){
            const doublec v = op(s, r);
            if (v == doublec(0.0)) continue;
            H(b, r ? (b | mask) : (b & ~mask)) += v;
        }
    }
}

// H += h4 (4x4 in the basis |s_i s_j>) on spins i != j
void addTwo(MatrixXcd& H, const Eigen::Matrix4cd& h4, int i, int j, int n){
    const int d = 1 << n, mi = 1 << (n - 1 - i), mj = 1 << (n - 1 - j);
    for (int b = 0; b < d; ++b){
        const int row = 2 * bitOf(b, i, n) + bitOf(b, j, n);
        const int base = b & ~mi & ~mj;
        for (int c = 0; c < 4; ++c){
            const doublec v = h4(row, c);
            if (v == doublec(0.0)) continue;
            H(b, base | ((c >> 1) ? mi : 0) | ((c & 1) ? mj : 0)) += v;
        }
    }
}

Eigen::Matrix4cd kron2(const Eigen::Matrix2cd& a, const Eigen::Matrix2cd& b){
    Eigen::Matrix4cd out;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) out.block<2, 2>(2 * i, 2 * j) = a(i, j) * b;
    return out;
}

// U <- (kron_k m[k]) U, spin 0 the first factor: one 2x2 per spin on the rows, O(n d^2) instead of a d x d product
void applyLocal(const std::vector<Eigen::Matrix2cd>& m, MatrixXcd& U, int c0 = 0, int c1 = -1){
    const int n = (int)m.size(), d = (int)U.rows();
    if (c1 < 0) c1 = (int)U.cols();
    for (int k = 0; k < n; ++k){
        const int mask = 1 << (n - 1 - k);
        const Eigen::Matrix2cd& a = m[k];
        for (int b = 0; b < d; ++b){
            if (b & mask) continue;
            const int b1 = b | mask;
            for (int c = c0; c < c1; ++c){
                const doublec u0 = U(b, c), u1 = U(b1, c);
                U(b, c) = a(0, 0) * u0 + a(0, 1) * u1;
                U(b1, c) = a(1, 0) * u0 + a(1, 1) * u1;
            }
        }
    }
}

// spin-1/2 operators J = sigma / 2
const Eigen::Matrix2cd& Jx(){ static Eigen::Matrix2cd m = (Eigen::Matrix2cd() << 0, 0.5, 0.5, 0).finished(); return m; }
const Eigen::Matrix2cd& Jy(){ static Eigen::Matrix2cd m = (Eigen::Matrix2cd() << 0, doublec(0, -0.5), doublec(0, 0.5), 0).finished(); return m; }
const Eigen::Matrix2cd& Jz(){ static Eigen::Matrix2cd m = (Eigen::Matrix2cd() << 0.5, 0, 0, -0.5).finished(); return m; }
Eigen::Matrix2cd along(const Eigen::Vector3d& v){ return v(0) * Jx() + v(1) * Jy() + v(2) * Jz(); }

Eigen::Matrix3d pairTensor(BathSpin* a, BathSpin* b){
    return calPointDipoleTensor(BathSpin_getXyz(a), BathSpin_getXyz(b), BathSpin_getGyro(a), BathSpin_getGyro(b)).real();
}

void checkFrames(const std::vector<Frame>& fr, int i, int j, double df, double Tzz){
    const double Wmax = std::max(std::hypot(fr[i].wx, fr[i].wy), std::hypot(fr[j].wx, fr[j].wy));
    if (Wmax > 0.0 && df < 3.0 * std::max(Wmax, std::fabs(Tzz))){
        char key[64], msg[256];
        snprintf(key, sizeof key, "frames %.3f", df);
        snprintf(msg, sizeof msg, "two different frames only %.3f MHz apart near a driven frame; their flip-flop is dropped "
                 "although it may be (Hartmann-Hahn) resonant", df / MHZ_TO_RADKHZ(1.0));
        warnOnce(key, msg);
    }
}

// Drive-assisted cross-line flip-flop: tone k tilts an off-line spin j by W_k / D_jk; in frame f_k the dropped
// flip-flop becomes dH = -(T_zz / 2) (W_k / D_jk) (cos p_k J1_k + sin p_k J2_k) Jn_j  (ops: per spin, the frame axes)
void addXff(MatrixXcd& H, BathArray* ba, const std::vector<Frame>& fr, const std::vector<std::array<Eigen::Matrix2cd, 3>>& ops,
            int n){
    for (int k = 0; k < n; ++k){
        if (fr[k].wx == 0.0 && fr[k].wy == 0.0) continue;
        const double Wk = std::hypot(fr[k].wx, fr[k].wy);
        for (int j = 0; j < n; ++j){
            const double Djk = BathSpin_getDetuning(ba->bath[j]) - fr[k].f;
            if (j == k || std::fabs(fr[j].f - fr[k].f) < kFrameEps || std::fabs(Djk) < 3.0 * Wk) continue;
            const double Tzz = pairTensor(ba->bath[k], ba->bath[j])(2, 2);
            addTwo(H, kron2(-0.5 * Tzz / Djk * (fr[k].wx * ops[k][0] + fr[k].wy * ops[k][1]), ops[j][2]), k, j, n);
        }
    }
}

// Per bath spin, the axes in qubit branch ms with the qubit coupling kept exactly: w_k z + ms A_k = W_k n_k
Axes branchAxes(BathArray* ba, Config* cnf, double ms){
    const double bz = fieldZ(cnf);
    const int n = BathArray_getNspin(ba);
    Axes ax;
    ax.n.assign(n, Eigen::Vector3d(0, 0, 1)); ax.e1.assign(n, Eigen::Vector3d(1, 0, 0)); ax.e2.assign(n, Eigen::Vector3d(0, 1, 0));
    ax.delta.assign(n, 0.0);
    for (int k = 0; k < n; ++k){
        BathSpin* bs = ba->bath[k];
        const MatrixXcd& A = bs->hypf[0];
        const double a0 = A(2, 0).real(), a1 = A(2, 1).real(), a2 = A(2, 2).real();
        const double vz = -BathSpin_getGyro(bs) * bz + BathSpin_getDetuning(bs) + ms * a2;
        const double ap2 = ms * ms * (a0 * a0 + a1 * a1);
        if (ap2 == 0.0){ ax.delta[k] = ms * a2; continue; }
        const double Wk = std::sqrt(vz * vz + ap2);
        ax.delta[k] = ms * a2 + (vz > 0 ? ap2 / (Wk + vz) : Wk - vz);   // no cancellation
        const Eigen::Vector3d nk = Eigen::Vector3d(ms * a0, ms * a1, vz) / Wk;
        Eigen::Vector3d x = Eigen::Vector3d(1, 0, 0) - nk(0) * nk;
        if (x.norm() < 1e-6) x = Eigen::Vector3d(0, 1, 0) - nk(1) * nk;
        ax.n[k] = nk;
        ax.e1[k] = x / x.norm();
        ax.e2[k] = nk.cross(ax.e1[k]);
    }
    return ax;
}

// One branch in its tilted frames
MatrixXcd tiltedHamiltonian(BathArray* ba, const std::vector<Frame>& fr, const Axes& ax, bool shifts){
    const int n = BathArray_getNspin(ba), d = 1 << n;
    std::vector<std::array<Eigen::Matrix2cd, 3>> ops(n);              // along e1, e2, n
    for (int k = 0; k < n; ++k) ops[k] = {along(ax.e1[k]), along(ax.e2[k]), along(ax.n[k])};
    MatrixXcd H = MatrixXcd::Zero(d, d);
    for (int k = 0; k < n; ++k){
        // the linear drive along x: only its part off n_k co-rotates about n_k
        const double s = std::sqrt(std::max(0.0, 1.0 - ax.n[k](0) * ax.n[k](0)));
        addOne(H, (fr[k].dz + ax.delta[k]) * ops[k][2] + s * (fr[k].wx * ops[k][0] + fr[k].wy * ops[k][1]), k, n);
    }
    if (shifts) addXff(H, ba, fr, ops, n);
    for (int i = 0; i < n; ++i){
        Eigen::Matrix3d Ri; Ri << ax.e1[i], ax.e2[i], ax.n[i];        // columns = local axes
        for (int j = i + 1; j < n; ++j){
            Eigen::Matrix3d Rj; Rj << ax.e1[j], ax.e2[j], ax.n[j];
            const Eigen::Matrix3d T = Ri.transpose() * pairTensor(ba->bath[i], ba->bath[j]) * Rj;
            Eigen::Matrix4cd h4 = T(2, 2) * kron2(ops[i][2], ops[j][2]);
            const double df = std::fabs(fr[i].f - fr[j].f);
            if (df < kFrameEps){
                h4 += 0.5 * (T(0, 0) + T(1, 1)) * (kron2(ops[i][0], ops[j][0]) + kron2(ops[i][1], ops[j][1]));
                h4 += 0.5 * (T(0, 1) - T(1, 0)) * (kron2(ops[i][0], ops[j][1]) - kron2(ops[i][1], ops[j][0]));
            }else checkFrames(fr, i, j, df, T(2, 2));
            addTwo(H, h4, i, j, n);
        }
    }
    return H;
}

// The A_zz-only model (qubit_transverse false): both branches in the same frames, which then cancel
void untiltedHamiltonians(BathArray* ba, Config* cnf, const std::vector<Frame>& fr, bool shifts, double msa, double msb,
                          MatrixXcd& Ha, MatrixXcd& Hb){
    const double bz = fieldZ(cnf);
    const int n = BathArray_getNspin(ba), d = 1 << n;
    std::vector<std::array<Eigen::Matrix2cd, 3>> ops(n, {Jx(), Jy(), Jz()});
    MatrixXcd Hbath = MatrixXcd::Zero(d, d), Hz = MatrixXcd::Zero(d, d), Hvv = MatrixXcd::Zero(d, d);
    for (int k = 0; k < n; ++k) addOne(Hbath, fr[k].dz * Jz() + fr[k].wx * Jx() + fr[k].wy * Jy(), k, n);
    if (shifts) addXff(Hbath, ba, fr, ops, n);
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j){
        const double Tzz = pairTensor(ba->bath[i], ba->bath[j])(2, 2);
        Eigen::Matrix4cd h4 = Tzz * kron2(Jz(), Jz());
        const double df = std::fabs(fr[i].f - fr[j].f);
        if (df < kFrameEps) h4 += -0.5 * Tzz * (kron2(Jx(), Jx()) + kron2(Jy(), Jy()));
        else checkFrames(fr, i, j, df, Tzz);
        addTwo(Hbath, h4, i, j, n);
    }
    for (int k = 0; k < n; ++k){
        BathSpin* bs = ba->bath[k];
        const MatrixXcd& A = bs->hypf[0];
        addOne(Hz, A(2, 2).real() * Jz(), k, n);
        if (shifts){
            const double wk = -BathSpin_getGyro(bs) * bz + BathSpin_getDetuning(bs);
            addOne(Hvv, (std::norm(A(2, 0)) + std::norm(A(2, 1))) / (2.0 * wk) * Jz(), k, n);   // Van Vleck
        }
    }
    Ha = Hbath + msa * Hz + (msa * msa) * Hvv;
    Hb = Hbath + msb * Hz + (msb * msb) * Hvv;
}

// R_q(t)^dag R_p(t) on every spin, R_m(t) = exp(-i omega_k t n_k^m . J_k)
std::vector<Eigen::Matrix2cd> frameChange(const std::vector<double>& omega, const Axes& p, const Axes& q, double t){
    const int n = (int)omega.size();
    std::vector<Eigen::Matrix2cd> m(n);
    for (int k = 0; k < n; ++k){
        const double th = 0.5 * omega[k] * t;
        const Eigen::Matrix2cd Rp = std::cos(th) * Eigen::Matrix2cd::Identity() - doublec(0, std::sin(th)) * (2.0 * along(p.n[k]));
        const Eigen::Matrix2cd Rq = std::cos(th) * Eigen::Matrix2cd::Identity() - doublec(0, std::sin(th)) * (2.0 * along(q.n[k]));
        m[k] = Rq.adjoint() * Rp;
    }
    return m;
}

// X = V exp(-i e tau) V^dagger (zgemm from ccexBlas::gemmMin on, as the cCCE Hahn step)
void propagator(const Eigen::VectorXd& e, const MatrixXcd& V, double tau, MatrixXcd& scratch, MatrixXcd& X){
    const Eigen::VectorXcd ph = (e * doublec(0.0, -tau)).array().exp();
#ifdef CCEX_BLAS_DIRECT
    if (V.rows() >= ccexBlas::gemmMin){ ccexBlas::hahnProduct(V, ph, scratch, X); return; }
#endif
    scratch.noalias() = V * ph.asDiagonal();
    X.noalias() = scratch * V.adjoint();
}

// C = A B
void multiply(const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C){
#ifdef CCEX_BLAS_DIRECT
    if (A.rows() >= ccexBlas::gemmMin){ ccexBlas::product(A, B, C); return; }
#endif
    C.noalias() = A * B;
}

// ---- eigensystems -------------------------------------------------------------------------------------------------

#ifdef CCEX_BLAS_DIRECT
extern "C" void dsyevd_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*,
                        const int*, int*);
#endif

// Hermitian H; a real symmetric one (|Im| < 1e-12, as pyCCEX tests it) by the real solver, about 3x faster: the m_s = 0
// branch and the untilted model usually are
void hermitianEigen(const MatrixXcd& H, CcexEigensolver& es, Eigen::VectorXd& w, MatrixXcd& V){
    if (H.imag().cwiseAbs().maxCoeff() >= 1e-12){
        es.compute(H); w = es.eigenvalues(); V = es.eigenvectors(); return;
    }
    Eigen::MatrixXd A = H.real();
    const int n = (int)A.rows();
#ifdef CCEX_BLAS_DIRECT
    if (n >= ccexBlas::eigMin){
        w.resize(n);
        int info = 0, lw = -1, liw = -1, iq = 0; double wq = 0;
        dsyevd_("V", "L", &n, A.data(), &n, w.data(), &wq, &lw, &iq, &liw, &info);
        std::vector<double> work(std::max(1, (int)wq)); std::vector<int> iwork(std::max(1, iq));
        lw = (int)work.size(); liw = (int)iwork.size();
        dsyevd_("V", "L", &n, A.data(), &n, w.data(), work.data(), &lw, iwork.data(), &liw, &info);
        if (info != 0){ es.compute(H); w = es.eigenvalues(); V = es.eigenvectors(); return; }
        V = A.cast<doublec>();
        return;
    }
#endif
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> rs(A);
    w = rs.eigenvalues(); V = rs.eigenvectors().cast<doublec>();
}

// ---- the mean field (configuration.s_range) ----------------------------------------------------------------------

double tzz(const double* x1, double g1, const double* x2, double g2){
    double a[3] = {x1[0], x1[1], x1[2]}, b[3] = {x2[0], x2[1], x2[2]};
    return calPointDipoleTensor(a, b, g1, g2)(2, 2).real();
}

// The field coefficients a source j (axis Pj, m along Pj) puts on spin i: the Hartree field of the pair Hamiltonian the
// cluster uses, T_zz m_j [Pjz Jz_i - (Pjx Jx_i + Pjy Jy_i)/2] in the same frame, T_zz m_j Pjz Jz_i in different frames;
// "z": T_zz m_j Jz_i
std::array<double, 3> mfCoupling(double T, bool same, const double* Pj, bool zmode){
    if (zmode) return {0.0, 0.0, T};
    const double s = same ? 1.0 : 0.0;
    return {-0.5 * T * s * Pj[0], -0.5 * T * s * Pj[1], T * Pj[2]};
}

// unit axis of a spin's own driven-frame Hamiltonian (wx, wy, dz); z when it has none
void frameAxis(const Frame& f, double* P){
    const double v[3] = {f.wx, f.wy, f.dz};
    const double nrm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (nrm > 0){ P[0] = v[0] / nrm; P[1] = v[1] / nrm; P[2] = v[2] / nrm; }
    else { P[0] = 0.0; P[1] = 0.0; P[2] = 1.0; }
}

// spin-1/2 eigenstate of n.J with eigenvalue m, basis [up, down]
Eigen::Vector2cd spinState(const double* nv, double m){
    const double th = std::acos(std::max(-1.0, std::min(1.0, nv[2]))), ph = std::atan2(nv[1], nv[0]);
    if (m > 0) return Eigen::Vector2cd(std::cos(th / 2), std::exp(doublec(0, ph)) * std::sin(th / 2));
    return Eigen::Vector2cd(-std::exp(doublec(0, -ph)) * std::sin(th / 2), std::cos(th / 2));
}

std::string sampleName(const std::string& name, int s){
    const size_t at = name.find("{s}");
    if (at == std::string::npos) return name + std::to_string(s);
    return name.substr(0, at) + std::to_string(s) + name.substr(at + 3);
}

}   // namespace

// ---------------------------------------------------------------------------------------------------------------------

void BathDrive_check(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls){
    const BathDrive* dr = pls->bdrive;
    if (!dr) return;
    const char* method = Config_getMethod(cnf);
    if (strcasecmp(method, "cce") && strcasecmp(method, "pcce")){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive runs with method cCCE (cluster.mode standard or pcce), not %s\n",
                         inputMethodName(method));
        ccexRefuse();
    }
    if (strcasecmp(Config_getQuantity(cnf), "coherence")){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive gives the coherence only (quantity coherence)\n"); ccexRefuse();
    }
    if (dr->mean_field != "pfull" && dr->mean_field != "projected" && dr->mean_field != "z" && dr->mean_field != "none"){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field must be pfull, projected, z or none, not \"%s\"\n",
                         dr->mean_field.c_str());
        ccexRefuse();
    }
    if (QubitArray_getNqubit(qa) != 1){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive takes one qubit\n"); ccexRefuse();
    }
    const float* b = Config_getBfield(cnf);
    if (!(std::isfinite(b[0]) && std::isfinite(b[1]) && std::isfinite(b[2])) ||
        std::hypot((double)b[0], (double)b[1]) > 1e-9 * std::max(std::fabs((double)b[2]), 1.0) || b[2] <= 0){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive needs bfield along +z (the frames are z rotations and the shift "
                         "signs assume a positive electron Zeeman frequency)\n");
        ccexRefuse();
    }
    if (pls->bnpulse > 0){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive and pulse.bath (DEER pulses) together are not implemented\n"); ccexRefuse();
    }
    // the qubit: instantaneous pi pulses of the default Ramsey / Hahn / CPMG sequence
    const int N = Pulse_getNpulse(pls);
    double** seq = Pulse_getSequence(pls);
    bool cpmg = Pulse_getPulseTime(pls) == 0.0 && Pulse_getPulseDetuningFactor(pls) == 1.0;
    for (int i = 0; cpmg && i <= N; ++i){
        const double start = (i == 0) ? 0.0 : double(2 * i - 1) / double(2 * N);
        const double end = (i == N) ? 1.0 : double(2 * (i + 1) - 1) / double(2 * N);
        if (std::fabs(seq[i][0] - start) > 1e-12 || std::fabs(seq[i][1] - end) > 1e-12) cpmg = false;
        if (i < N){
            const char axis = toupper(Pulse_getPulseAxes(pls)[i]);
            if ((axis != 'X' && axis != 'Y') || std::fabs(std::fabs(Pulse_getPulseAngles(pls)[i]) - 180.0) > 1e-9) cpmg = false;
        }
    }
    if (!cpmg){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive takes the default sequence of pulse.count (Ramsey, Hahn echo, CPMG): "
                         "instantaneous X or Y pi pulses at (2c-1)/2N of the time, no duration, detuning_factor 1\n");
        ccexRefuse();
    }
    if (!(dr->line_tolerance > 0.0) || !std::isfinite(dr->line_tolerance)){
        ccexRefusalPrint(stderr, "Error: pulse.bath_drive.line_tolerance must be finite and > 0\n"); ccexRefuse();
    }
    const std::vector<Component> comps = components(*dr);
    for (size_t i = 0; i < comps.size(); ++i) for (size_t j = i + 1; j < comps.size(); ++j)
        if (std::fabs(comps[i].carrier - comps[j].carrier) < dr->line_tolerance){
            ccexRefusalPrint(stderr, "Error: pulse.bath_drive: two drive components share a carrier (%.6f MHz); combine them "
                             "into one tone\n", comps[i].carrier / MHZ_TO_RADKHZ(1.0));
            ccexRefuse();
        }
    const std::vector<Frame> fr = frameAssignment(ba, cnf, *dr, comps, true);
    if (rank == 0){
        printf("\n    >> Continuous bath drive (pulse.bath_drive)\n");
        printf("      shifts           : %s\n", dr->shifts ? "true" : "false");
        printf("      qubit_transverse : %s\n", dr->qubit_transverse ? "true" : "false");
        if (Config_getNstate(cnf) > 0)
            printf("      mean_field       : %s (sampled states: one pass, every cluster averaged over them)\n",
                   dr->mean_field.c_str());
        else if (!dr->shell_file.empty() || !dr->shell_state_file.empty())
            printf("      Note: mean_field_shell and mean_field_shell_state_file are not used without configuration.s_range\n");
        for (size_t it = 0; it < dr->tones.size(); ++it){
            const BathDriveTone& t = dr->tones[it];
            int driven = 0;
            for (int k = 0; k < (int)fr.size(); ++k)
                if (fr[k].wx != 0.0 || fr[k].wy != 0.0)
                    for (const Component& c : comps) if ((int)it == c.tone && std::fabs(fr[k].f - c.carrier) < kFrameEps){ ++driven; break; }
            printf("      tone %zu : rabi %.6f MHz, offset %.6f MHz, phase %.3f deg, lines [", it, t.rabi / MHZ_TO_RADKHZ(1.0),
                   t.offset / MHZ_TO_RADKHZ(1.0), t.phase * 180.0 / M_PI);
            for (size_t l = 0; l < t.lines.size(); ++l) printf("%s%.6f", l ? ", " : "", t.lines[l] / MHZ_TO_RADKHZ(1.0));
            printf("] MHz, %d of %d bath spins driven\n", driven, (int)fr.size());
        }
    }
}


void BathDrive_prepareMeanField(BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls){
    BathDrive& dr = *pls->bdrive;
    const int last = Config_getNstate(cnf), first = std::max(1, cnf->state_first);
    const int N = BathArray_getNspin(ba), ns = last - first + 1;
    const bool zmode = dr.mean_field == "z";
    dr.ns = ns; dr.nbath = N;
    // the samples: the bath states of s = first .. last as the run would set them (bath_state_file or the seed); the
    // states below first are drawn, as the state loop draws them. m is read along each spin's axis P
    dr.m.assign((size_t)ns * N, 0.0);
    for (int s = 1; s <= last; ++s){
        if (s < first){ drawSkippedState(ba, dfa, cnf); continue; }
        setBathStates(ba, cnf, s);
        for (int i = 0; i < N; ++i){
            const double m = BathArray_getBath_i_state(ba, i);
            if (std::fabs(std::fabs(m) - 0.5) > 1e-12){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive: the state of bath spin %d in sample %d is %g, not +-1/2\n", i, s, m);
                ccexRefuse();
            }
            dr.m[(size_t)(s - first) * N + i] = m;
        }
    }
    const std::vector<Component> comps = components(dr);
    const std::vector<Frame> fr = frameAssignment(ba, cnf, dr, comps, false);
    dr.P.assign((size_t)3 * N, 0.0);
    for (int i = 0; i < N; ++i){
        if (zmode) dr.P[3 * i + 2] = 1.0;
        else frameAxis(fr[i], &dr.P[3 * i]);
    }

    // the shell: spins that only give a field, "x y z axis [name]" (A), a lone number (a count) or # skipped
    std::vector<std::array<double, 3>> xyzS;
    std::vector<SpinData> spS;
    if (!dr.shell_file.empty()){
        std::ifstream in(dr.shell_file);
        if (!in){
            ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell: cannot read %s\n", dr.shell_file.c_str());
            ccexRefuse();
        }
        std::string line;
        int ln = 0;
        while (std::getline(in, line)){
            ++ln;
            std::istringstream is(line);
            std::vector<std::string> tok;
            for (std::string w; is >> w;) tok.push_back(w);
            if (tok.empty() || tok[0][0] == '#' || tok.size() == 1) continue;
            if (tok.size() != 4 && tok.size() != 5){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell line %d: expected \"x y z axis [name]\"\n", ln);
                ccexRefuse();
            }
            const std::string name = tok.size() == 5 ? tok[4] : dr.shell_species;
            const int ax = (int)std::stod(tok[3]);
            int ib = -1;
            for (int i = 0; i < N && ib < 0; ++i) if (!strcasecmp(BathArray_getBath_i_name(ba, i), name.c_str())) ib = i;
            if (ib < 0){
                ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell line %d: species %s is not in the bath (its spin "
                                 "and gyro are taken from there)\n", ln, name.c_str());
                ccexRefuse();
            }
            const int idf = DefectArray_findDefectIndex(dfa, (char*)name.c_str());
            double det = 0.0;
            if (idf < 0){
                if (ax != 0){
                    ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell line %d: axis %d but no Defect table for %s\n",
                                     ln, ax, name.c_str());
                    ccexRefuse();
                }
            }else{
                if (ax < 1 || ax > DefectArray_getDefect_idf_navaax(dfa, idf)){
                    ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell line %d: axis %d is not in the %s table\n",
                                     ln, ax, name.c_str());
                    ccexRefuse();
                }
                det = MHZ_TO_RADKHZ(DefectArray_getDefect_idf_iax_detuning(dfa, idf, ax));
            }
            xyzS.push_back({std::stod(tok[0]), std::stod(tok[1]), std::stod(tok[2])});
            spS.push_back({BathArray_getBath_i_gyro(ba, ib), det, (double)BathArray_getBath_i_spin(ba, ib), name});
        }
    }
    const int NS = (int)xyzS.size();
    const std::vector<Frame> frS = frameAssignment(spS, cnf, dr, comps, true);   // refuses a spin that is not spin-1/2 electron-like
    std::vector<double> PS((size_t)3 * NS, 0.0);
    for (int j = 0; j < NS; ++j){
        if (zmode) PS[3 * j + 2] = 1.0;
        else frameAxis(frS[j], &PS[3 * j]);
    }
    // the shell's samples: mean_field_shell_state_file per sample, else drawn from seed + 7919 (every rank the same)
    std::vector<double> mS((size_t)ns * NS, 0.0);
    if (NS){
        std::mt19937_64 rng((unsigned long long)Config_getSeed(cnf) + 7919ULL);
        for (int s = 1; s <= last; ++s){
            std::vector<double> row(NS);
            if (!dr.shell_state_file.empty()){
                if (s < first) continue;
                const std::string fn = sampleName(dr.shell_state_file, s);
                std::ifstream in(fn);
                std::string skip;
                if (!in || !std::getline(in, skip)){
                    ccexRefusalPrint(stderr, "Error: pulse.bath_drive.mean_field_shell_state_file: cannot read %s\n", fn.c_str());
                    ccexRefuse();
                }
                for (int j = 0; j < NS; ++j)
                    if (!(in >> row[j]) || std::fabs(std::fabs(row[j]) - 0.5) > 1e-12){
                        ccexRefusalPrint(stderr, "Error: %s: line %d of the shell must hold +-1/2\n", fn.c_str(), j + 2);
                        ccexRefuse();
                    }
            }else for (int j = 0; j < NS; ++j) row[j] = (rng() & 1ULL) ? 0.5 : -0.5;
            if (s >= first) std::copy(row.begin(), row.end(), mS.begin() + (size_t)(s - first) * NS);
        }
    }

    // h = the field of every other bath spin and of the shell, per sample
    dr.h.assign((size_t)ns * N * 3, 0.0);
    for (int i = 0; i < N; ++i){
        double* xi = BathArray_getBath_i_xyz(ba, i);
        const double gi = BathArray_getBath_i_gyro(ba, i);
        auto add = [&](const std::array<double, 3>& g, const double* mj, int stride){
            for (int s = 0; s < ns; ++s){
                const double m = mj[(size_t)s * stride];
                double* hs = &dr.h[((size_t)s * N + i) * 3];
                hs[0] += g[0] * m; hs[1] += g[1] * m; hs[2] += g[2] * m;
            }
        };
        for (int j = 0; j < N; ++j){
            if (j == i) continue;
            const double T = tzz(xi, gi, BathArray_getBath_i_xyz(ba, j), BathArray_getBath_i_gyro(ba, j));
            add(mfCoupling(T, std::fabs(fr[i].f - fr[j].f) < kFrameEps, &dr.P[3 * j], zmode), &dr.m[j], N);
        }
        for (int j = 0; j < NS; ++j){
            const double T = tzz(xi, gi, xyzS[j].data(), spS[j].gyro);
            add(mfCoupling(T, std::fabs(fr[i].f - frS[j].f) < kFrameEps, &PS[3 * j], zmode), &mS[j], NS);
        }
    }
    if (rank == 0)
        printf("      mean field       : %d samples (s = %d .. %d), %d bath spins, %d shell spins%s\n\n", ns, first, last, N, NS,
               dr.shell_state_file.empty() && NS ? " (shell states from seed + 7919)" : "");
}

MatrixXcd* calCoherenceCceDriven(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, const int* idx){
    const BathDrive& dr = *pls->bdrive;
    const std::vector<Component> comps = components(dr);
    const std::vector<Frame> fr = frameAssignment(ba, cnf, dr, comps, false);
    const double msa = findZbasisSubLevel(qa->psia), msb = findZbasisSubLevel(qa->psib);
    const int n = BathArray_getNspin(ba), d = 1 << n;

    MatrixXcd Ha, Hb;
    Axes axa, axb;
    std::vector<double> omega;
    const bool tilted = dr.qubit_transverse;
    if (tilted){
        axa = branchAxes(ba, cnf, msa);
        axb = branchAxes(ba, cnf, msb);
        Ha = tiltedHamiltonian(ba, fr, axa, dr.shifts);
        Hb = tiltedHamiltonian(ba, fr, axb, dr.shifts);
        omega.resize(n);
        for (int k = 0; k < n; ++k) omega[k] = -BathSpin_getGyro(ba->bath[k]) * fieldZ(cnf) + fr[k].f;
    }else untiltedHamiltonians(ba, cnf, fr, dr.shifts, msa, msb, Ha, Hb);

    const int nstep = Config_getNstep(cnf);
    const double deltat = Config_getDeltat(cnf);
    const int N = Pulse_getNpulse(pls);

    if (dr.ns > 0){
        // Sampled states: per sample the cluster starts in the product state along the axes P (m = +-1/2), the field of
        // every spin outside the cluster is static (the mean field less the cluster's own spins), and the coherence is
        // the average over the samples (internal averaging)
        const int ns = dr.ns, NB = dr.nbath;
        const bool zmode = dr.mean_field == "z", projected = dr.mean_field == "projected", none = dr.mean_field == "none";
        std::vector<std::array<double, 3>> g((size_t)n * n, {0.0, 0.0, 0.0});
        for (int k = 0; k < n; ++k) for (int l = 0; l < n; ++l) if (k != l){
            const double T = pairTensor(ba->bath[k], ba->bath[l])(2, 2);
            g[(size_t)k * n + l] = mfCoupling(T, std::fabs(fr[k].f - fr[l].f) < kFrameEps, &dr.P[3 * idx[l]], zmode);
        }
        std::vector<std::array<Eigen::Matrix2cd, 3>> opsA(n), opsB(n);   // the field's axes: the branch's tilted axes, or x y z
        for (int k = 0; k < n; ++k){
            if (tilted){
                opsA[k] = {along(axa.e1[k]), along(axa.e2[k]), along(axa.n[k])};
                opsB[k] = {along(axb.e1[k]), along(axb.e2[k]), along(axb.n[k])};
            }else opsA[k] = opsB[k] = {Jx(), Jy(), Jz()};
        }
        Eigen::VectorXd t(nstep);
        for (int i = 0; i < nstep; ++i) t(i) = i * deltat;
        Eigen::VectorXcd L = Eigen::VectorXcd::Zero(nstep);
        MatrixXcd Hsa, Hsb, VaH, VbH, va, vb, tmp, Va, Vb;
        Eigen::VectorXd ea, eb;
        CcexEigensolver sa, sb;
        for (int s = 0; s < ns; ++s){
            Hsa = Ha; Hsb = Hb;
            for (int k = 0; k < n; ++k){
                double h[3];
                for (int c = 0; c < 3; ++c){
                    h[c] = dr.h[((size_t)s * NB + idx[k]) * 3 + c];
                    for (int l = 0; l < n; ++l) if (l != k) h[c] -= g[(size_t)k * n + l][c] * dr.m[(size_t)s * NB + idx[l]];
                }
                const double* Pk = &dr.P[3 * idx[k]];
                if (projected){
                    const double dot = h[0] * Pk[0] + h[1] * Pk[1] + h[2] * Pk[2];
                    for (int c = 0; c < 3; ++c) h[c] = dot * Pk[c];
                }else if (none) h[0] = h[1] = h[2] = 0.0;
                addOne(Hsa, h[0] * opsA[k][0] + h[1] * opsA[k][1] + h[2] * opsA[k][2], k, n);
                addOne(Hsb, h[0] * opsB[k][0] + h[1] * opsB[k][1] + h[2] * opsB[k][2], k, n);
            }
            hermitianEigen(Hsa, sa, ea, Va);
            hermitianEigen(Hsb, sb, eb, Vb);
            VaH = Va.adjoint(); VbH = Vb.adjoint();
            // the initial product state, the same in every time column
            Eigen::VectorXcd psi = Eigen::VectorXcd::Ones(1);
            for (int k = 0; k < n; ++k){
                const Eigen::Vector2cd st = spinState(&dr.P[3 * idx[k]], dr.m[(size_t)s * NB + idx[k]]);
                Eigen::VectorXcd next(psi.size() * 2);
                for (int r = 0; r < psi.size(); ++r){ next(2 * r) = psi(r) * st(0); next(2 * r + 1) = psi(r) * st(1); }
                psi.swap(next);
            }
            va = psi.replicate(1, nstep); vb = va;
            // v <- V exp(-i e tau) V^dagger v, column i with its own tau_i
            auto apply = [&](bool branchA, const Eigen::VectorXd& tau, MatrixXcd& v){
                const Eigen::VectorXd& e = branchA ? ea : eb;
                multiply(branchA ? VaH : VbH, v, tmp);
                for (int i = 0; i < nstep; ++i) tmp.col(i).array() *= (e * doublec(0.0, -tau(i))).array().exp();
                multiply(branchA ? Va : Vb, tmp, v);
            };
            if (N == 0){ apply(true, t, va); apply(false, t, vb); }
            else{
                const Eigen::VectorXd h1 = t / (2.0 * N);
                apply(true, h1, va); apply(false, h1, vb);
                for (int k = 1; k <= N; ++k){
                    const bool odd = k % 2 == 1;
                    const Eigen::VectorXd tau = (k == N) ? h1 : Eigen::VectorXd(2.0 * h1);
                    if (tilted){
                        const Axes& oldA = odd ? axa : axb;
                        const Axes& newA = odd ? axb : axa;
                        for (int i = 0; i < nstep; ++i){
                            const double sk = t(i) / (2 * N) + (k - 1) * t(i) / N;
                            applyLocal(frameChange(omega, oldA, newA, sk), va, i, i + 1);
                            applyLocal(frameChange(omega, newA, oldA, sk), vb, i, i + 1);
                        }
                    }
                    apply(!odd, tau, va);                         // path a: beta after an odd number of pi pulses
                    apply(odd, tau, vb);
                }
            }
            if (tilted){
                const Axes& endA = (N % 2 == 0) ? axa : axb;
                const Axes& endB = (N % 2 == 0) ? axb : axa;
                for (int i = 0; i < nstep; ++i) applyLocal(frameChange(omega, endA, endB, t(i)), va, i, i + 1);
            }
            L += (vb.conjugate().array() * va.array()).colwise().sum().transpose().matrix();
        }
        L /= doublec(ns, 0.0);
        MatrixXcd* result = new MatrixXcd[nstep];
        for (int i = 0; i < nstep; ++i){
            const doublec v = (N >= 2 && N % 2 == 0) ? std::conj(L(i)) : L(i);
            result[i] = MatrixXcd::Constant(1, 1, v);
        }
        return result;
    }

    CcexEigensolver sa, sb;                                      // zheevd from ccexBlas::eigMin on
    sa.compute(Ha); sb.compute(Hb);
    const Eigen::VectorXd& ea = sa.eigenvalues(); const MatrixXcd& Va = sa.eigenvectors();
    const Eigen::VectorXd& eb = sb.eigenvalues(); const MatrixXcd& Vb = sb.eigenvectors();

    MatrixXcd* result = new MatrixXcd[nstep];
    MatrixXcd Ua, Ub, halfA, halfB, fullA, fullB, scratch, next;
    for (int i = 0; i < nstep; ++i){
        const double t = i * deltat;
        if (N == 0){
            propagator(ea, Va, t, scratch, Ua);
            propagator(eb, Vb, t, scratch, Ub);
        }else{
            propagator(ea, Va, t / (2 * N), scratch, halfA);
            propagator(eb, Vb, t / (2 * N), scratch, halfB);
            if (N >= 2){ multiply(halfA, halfA, fullA); multiply(halfB, halfB, fullB); }
            Ua = halfA; Ub = halfB;                              // the first interval: a in alpha, b in beta
            for (int k = 1; k <= N; ++k){                        // after the k-th pi the branches swap
                const bool odd = k % 2 == 1, last = k == N;
                const MatrixXcd& na = odd ? (last ? halfB : fullB) : (last ? halfA : fullA);
                const MatrixXcd& nb = odd ? (last ? halfA : fullA) : (last ? halfB : fullB);
                if (tilted){                                     // frame change at s_k = t/2N + (k-1) t/N
                    const double sk = t / (2 * N) + (k - 1) * t / N;
                    const Axes& oldA = odd ? axa : axb;          // path a: from its old branch into the new one
                    const Axes& newA = odd ? axb : axa;
                    applyLocal(frameChange(omega, oldA, newA, sk), Ua);
                    applyLocal(frameChange(omega, newA, oldA, sk), Ub);
                }
                multiply(na, Ua, next); Ua.swap(next);
                multiply(nb, Ub, next); Ub.swap(next);
            }
        }
        if (tilted){                                             // R_endB(t)^dag R_endA(t) before the overlap
            const Axes& endA = (N % 2 == 0) ? axa : axb;
            const Axes& endB = (N % 2 == 0) ? axb : axa;
            applyLocal(frameChange(omega, endA, endB, t), Ua);
        }
        doublec L = (Ua.array() * Ub.conjugate().array()).sum() / doublec(d, 0.0);
        // the kernel's convention: Tr[Ua Ub^dag]/d is 2<beta|rho|alpha> for odd N and 2<alpha|rho|beta> for even N,
        // so even N >= 2 is conjugated (Ramsey: 2<alpha|rho|beta>); Re L is the same either way
        if (N >= 2 && N % 2 == 0) L = std::conj(L);
        result[i] = MatrixXcd::Constant(1, 1, L);
    }
    return result;
}
