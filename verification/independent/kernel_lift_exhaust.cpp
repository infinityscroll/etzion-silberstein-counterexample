#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

using u16 = std::uint16_t;
using u64 = std::uint64_t;

static int rank4(u16 a) {
    std::array<unsigned, 4> rows{};
    for (int i = 0; i < 4; ++i) rows[i] = (a >> (4 * i)) & 15u;
    int r = 0;
    for (int c = 0; c < 4; ++c) {
        int p = r;
        while (p < 4 && ((rows[p] >> c) & 1u) == 0) ++p;
        if (p == 4) continue;
        std::swap(rows[p], rows[r]);
        for (int i = 0; i < 4; ++i)
            if (i != r && ((rows[i] >> c) & 1u)) rows[i] ^= rows[r];
        ++r;
    }
    return r;
}

static u16 transpose4(u16 a) {
    u16 z = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            z |= ((a >> (4 * i + j)) & 1u) << (4 * j + i);
    return z;
}

static u16 multiply4(u16 a, u16 b) {
    u16 z = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned ar = (a >> (4 * i)) & 15u;
        unsigned zr = 0;
        for (int j = 0; j < 4; ++j) {
            unsigned bc = 0;
            for (int k = 0; k < 4; ++k) bc |= ((b >> (4 * k + j)) & 1u) << k;
            zr |= (std::popcount(ar & bc) & 1u) << j;
        }
        z |= zr << (4 * i);
    }
    return z;
}

static std::optional<u16> inverse4(u16 a) {
    std::array<unsigned, 4> r{};
    for (int i = 0; i < 4; ++i) r[i] = ((a >> (4 * i)) & 15u) | (1u << (4 + i));
    for (int c = 0; c < 4; ++c) {
        int p = c;
        while (p < 4 && ((r[p] >> c) & 1u) == 0) ++p;
        if (p == 4) return std::nullopt;
        std::swap(r[p], r[c]);
        for (int i = 0; i < 4; ++i) if (i != c && ((r[i] >> c) & 1u)) r[i] ^= r[c];
    }
    u16 z = 0;
    for (int i = 0; i < 4; ++i) z |= ((r[i] >> 4) & 15u) << (4 * i);
    return z;
}

static std::vector<u16> dual_basis(std::vector<u16> equations) {
    int r = 0;
    std::vector<int> pivots;
    for (int c = 0; c < 16 && r < static_cast<int>(equations.size()); ++c) {
        int p = r;
        while (p < static_cast<int>(equations.size()) && ((equations[p] >> c) & 1u) == 0) ++p;
        if (p == static_cast<int>(equations.size())) continue;
        std::swap(equations[p], equations[r]);
        for (int i = 0; i < static_cast<int>(equations.size()); ++i)
            if (i != r && ((equations[i] >> c) & 1u)) equations[i] ^= equations[r];
        pivots.push_back(c);
        ++r;
    }
    std::array<bool, 16> pivot{};
    for (int p : pivots) pivot[p] = true;
    std::vector<u16> basis;
    for (int f = 0; f < 16; ++f) if (!pivot[f]) {
        u16 x = u16(1u << f);
        for (int i = 0; i < r; ++i)
            if ((equations[i] >> f) & 1u) x |= u16(1u << pivots[i]);
        basis.push_back(x);
    }
    return basis;
}

static int binary_rank(std::vector<u16> rows, int width) {
    int r = 0;
    for (int c = 0; c < width; ++c) {
        int p = r;
        while (p < static_cast<int>(rows.size()) && ((rows[p] >> c) & 1u) == 0) ++p;
        if (p == static_cast<int>(rows.size())) continue;
        std::swap(rows[p], rows[r]);
        for (int i = 0; i < static_cast<int>(rows.size()); ++i)
            if (i != r && ((rows[i] >> c) & 1u)) rows[i] ^= rows[r];
        ++r;
    }
    return r;
}

static u16 encode(unsigned message, const std::vector<u16>& basis) {
    u16 a = 0;
    for (unsigned i = 0; i < basis.size(); ++i)
        if ((message >> i) & 1u) a ^= basis[i];
    return a;
}

static std::vector<u16> code_span(const std::vector<u16>& basis) {
    std::vector<u16> out(1u << basis.size());
    for (unsigned x = 0; x < out.size(); ++x) out[x] = encode(x, basis);
    return out;
}

struct RankTwoWord {
    u16 message;
    std::uint8_t annihilator_a;
    std::uint8_t annihilator_b;
};

static RankTwoWord describe_rank_two(u16 message, u16 matrix) {
    std::vector<unsigned> annihilator;
    for (unsigned v = 0; v < 16; ++v) {
        bool ok = true;
        for (int i = 0; i < 4; ++i) {
            unsigned row = (matrix >> (4 * i)) & 15u;
            if (std::popcount(v & row) & 1u) { ok = false; break; }
        }
        if (ok) annihilator.push_back(v);
    }
    if (annihilator.size() != 4) {
        std::cerr << "internal error: rank-two annihilator has size " << annihilator.size() << "\n";
        std::exit(3);
    }
    return RankTwoWord{message, static_cast<std::uint8_t>(annihilator[1]),
                       static_cast<std::uint8_t>(annihilator[2])};
}

struct Kernel {
    u16 h1, h2;
    std::uint8_t p, q;
};

static std::vector<Kernel> enumerate_kernels() {
    std::vector<Kernel> out;
    out.reserve(2794155);
    for (int p = 0; p < 11; ++p) for (int q = p + 1; q < 12; ++q) {
        const int exponent = (q - p - 1) + 2 * (11 - q);
        const unsigned limit = 1u << exponent;
        for (unsigned code = 0; code < limit; ++code) {
            unsigned t = code;
            u16 h1 = u16(1u << p), h2 = u16(1u << q);
            for (int j = p + 1; j < q; ++j) {
                if (t & 1u) h1 |= u16(1u << j);
                t >>= 1;
            }
            for (int j = q + 1; j < 12; ++j) {
                if (t & 1u) h1 |= u16(1u << j);
                if (t & 2u) h2 |= u16(1u << j);
                t >>= 2;
            }
            out.push_back(Kernel{h1, h2, static_cast<std::uint8_t>(p), static_cast<std::uint8_t>(q)});
        }
    }
    return out;
}

static u64 kernel_index_from_rows(u16 a, u16 b) {
    if (!a || !b || a == b) {
        std::cerr << "dependent annihilator rows\n";
        std::exit(3);
    }
    int p = std::countr_zero(u16(a | b));
    if (((a >> p) & 1u) == 0) std::swap(a, b);
    if ((b >> p) & 1u) b ^= a;
    int q = std::countr_zero(b);
    if ((a >> q) & 1u) a ^= b;
    u64 offset = 0;
    for (int pp = 0; pp < 11; ++pp) for (int qq = pp + 1; qq < 12; ++qq) {
        int exponent = (qq - pp - 1) + 2 * (11 - qq);
        if (pp == p && qq == q) {
            unsigned code = 0, shift = 0;
            for (int j = p + 1; j < q; ++j, ++shift)
                code |= ((a >> j) & 1u) << shift;
            for (int j = q + 1; j < 12; ++j, shift += 2) {
                code |= ((a >> j) & 1u) << shift;
                code |= ((b >> j) & 1u) << (shift + 1);
            }
            return offset + code;
        }
        offset += 1ull << exponent;
    }
    std::cerr << "RREF index failure\n";
    std::exit(3);
}

struct MatrixAutomorphism { u16 p, q; };
struct FunctionalAction { std::array<u16, 4096> image{}; };

static std::vector<MatrixAutomorphism> spread_automorphisms(
        const std::vector<u16>& generators) {
    auto spread = code_span(generators);
    std::array<bool, 65536> in_spread{};
    for (u16 a : spread) in_spread[a] = true;
    std::vector<MatrixAutomorphism> autos;
    for (unsigned q0 = 0; q0 < 65536; ++q0) {
        u16 q = static_cast<u16>(q0);
        auto qi = inverse4(q);
        if (!qi) continue;
        for (u16 b : spread) if (b) {
            u16 p = multiply4(b, *qi); // Force P I Q = b.
            bool ok = true;
            for (u16 g : generators) {
                if (!in_spread[multiply4(multiply4(p, g), q)]) { ok = false; break; }
            }
            if (ok) autos.push_back({p, q});
        }
    }
    std::unordered_set<u64> keys;
    for (auto a : autos) keys.insert(u64(a.p) | (u64(a.q) << 16));
    if (keys.size() != autos.size()) {
        std::cerr << "duplicate automorphisms\n";
        std::exit(3);
    }
    for (auto a : autos) for (auto b : autos) {
        // Composition A -> Pa(Pb A Qb)Qa.
        u16 p = multiply4(a.p, b.p), q = multiply4(b.q, a.q);
        if (!keys.contains(u64(p) | (u64(q) << 16))) {
            std::cerr << "automorphism set is not closed\n";
            std::exit(3);
        }
    }
    return autos;
}

static std::vector<FunctionalAction> functional_actions(
        const std::vector<MatrixAutomorphism>& autos, const std::vector<u16>& dual) {
    std::array<u16, 65536> matrix_message{};
    std::array<bool, 65536> in_code{};
    for (unsigned x = 0; x < 4096; ++x) {
        u16 a = encode(x, dual);
        matrix_message[a] = x;
        in_code[a] = true;
    }
    std::vector<FunctionalAction> actions;
    actions.reserve(autos.size());
    for (auto automorphism : autos) {
        u16 lp = transpose4(*inverse4(automorphism.p));
        u16 rq = transpose4(*inverse4(automorphism.q));
        std::array<u16, 4096> msg_image{}, msg_inverse{};
        for (unsigned x = 0; x < 4096; ++x) {
            u16 a = multiply4(multiply4(lp, encode(x, dual)), rq);
            if (!in_code[a]) {
                std::cerr << "dual automorphism failure\n";
                std::exit(3);
            }
            msg_image[x] = matrix_message[a];
            msg_inverse[msg_image[x]] = x;
        }
        FunctionalAction f;
        for (unsigned h = 0; h < 4096; ++h) {
            u16 hp = 0;
            for (int j = 0; j < 12; ++j)
                hp |= (std::popcount(u16(h & msg_inverse[1u << j])) & 1u) << j;
            f.image[h] = hp;
        }
        actions.push_back(f);
    }
    return actions;
}

struct OrbitRepresentative { u64 index; unsigned mass; };

static std::vector<OrbitRepresentative> kernel_orbits(
        const std::vector<Kernel>& kernels, const std::vector<FunctionalAction>& actions) {
    std::vector<bool> visited(kernels.size());
    std::vector<OrbitRepresentative> reps;
    u64 mass = 0;
    std::vector<u64> orbit;
    orbit.reserve(actions.size());
    for (u64 i = 0; i < kernels.size(); ++i) if (!visited[i]) {
        orbit.clear();
        for (const auto& action : actions)
            orbit.push_back(kernel_index_from_rows(action.image[kernels[i].h1],
                                                   action.image[kernels[i].h2]));
        std::sort(orbit.begin(), orbit.end());
        orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
        if (!std::binary_search(orbit.begin(), orbit.end(), i)) {
            std::cerr << "orbit omits representative\n";
            std::exit(3);
        }
        for (u64 j : orbit) {
            if (visited[j]) {
                std::cerr << "overlapping orbit partition\n";
                std::exit(3);
            }
            visited[j] = true;
        }
        reps.push_back({i, static_cast<unsigned>(orbit.size())});
        mass += orbit.size();
    }
    if (mass != kernels.size()) {
        std::cerr << "orbit mass failure\n";
        std::exit(3);
    }
    return reps;
}

static u16 compress_free(u16 x, int p, int q) {
    u16 z = 0;
    int k = 0;
    for (int j = 0; j < 12; ++j) {
        if (j == p || j == q) continue;
        z |= ((x >> j) & 1u) << k++;
    }
    return z;
}

struct Constraint { u64 e1, e2; };

// Affine equations over 48 Boolean unknowns: the four coordinates of the image
// of each of the 12 fixed basis vectors of U.  Any map on K extends to U, so
// this deliberately avoids choosing a basis of K.  Row p has least pivot p.
static constexpr int NVARIABLES = 48;
struct AffineSystem {
    std::array<u64, NVARIABLES> row{};
    std::array<std::uint8_t, NVARIABLES> rhs{};

    bool add(u64 coeff, bool value) {
        while (coeff) {
            int p = std::countr_zero(coeff);
            if (row[p]) {
                coeff ^= row[p];
                value ^= rhs[p];
            } else {
                row[p] = coeff;
                rhs[p] = value;
                return true;
            }
        }
        return !value;
    }

    struct Evaluation { bool fixed; bool value; u64 remainder; bool offset; };
    Evaluation evaluate(u64 coeff) const {
        bool value = false;
        u64 remainder = 0;
        while (coeff) {
            int p = std::countr_zero(coeff);
            if (row[p]) {
                coeff ^= row[p];
                value ^= rhs[p];
            } else {
                remainder |= 1ull << p;
                coeff ^= 1ull << p;
            }
        }
        return Evaluation{remainder == 0, value, remainder, value};
    }
};

struct SolveStats { u64 nodes = 0, conflicts = 0, forced = 0; };

struct ActiveSet {
    // A kernel has at most 525 rank-two words, hence nine 64-bit words suffice.
    std::array<u64, 9> bits{};
    void clear(unsigned i) { bits[i / 64] &= ~(1ull << (i % 64)); }
    bool empty() const {
        for (u64 x : bits) if (x) return false;
        return true;
    }
};

static bool solve_active(const std::vector<Constraint>& constraints,
                         AffineSystem system, SolveStats& stats, ActiveSet active) {
    ++stats.nodes;
    for (;;) {
        int branch = -1;
        int best_score = std::numeric_limits<int>::max();
        bool restarted = false;
        ActiveSet remaining = active;
        for (unsigned word = 0; word < active.bits.size() && !restarted; ++word) {
          u64 pending = active.bits[word];
          while (pending) {
            unsigned bit = std::countr_zero(pending);
            pending &= pending - 1;
            const unsigned i = 64 * word + bit;
            const auto a = system.evaluate(constraints[i].e1);
            const auto b = system.evaluate(constraints[i].e2);
            if ((a.fixed && a.value) || (b.fixed && b.value)) {
                remaining.clear(i);
                continue;
            }
            if (a.fixed && b.fixed) {
                ++stats.conflicts;
                return false;
            }
            if (a.fixed && !a.value) {
                if (!system.add(constraints[i].e2, true)) {
                    ++stats.conflicts;
                    return false;
                }
                ++stats.forced;
                remaining.clear(i);
                active = remaining;
                restarted = true;
                break;
            }
            if (b.fixed && !b.value) {
                if (!system.add(constraints[i].e1, true)) {
                    ++stats.conflicts;
                    return false;
                }
                ++stats.forced;
                remaining.clear(i);
                active = remaining;
                restarted = true;
                break;
            }
            // Equal reduced linear parts mean the two forms are either equal
            // (a unit equation) or complementary (the clause is automatic).
            if (a.remainder == b.remainder) {
                if (a.offset != b.offset) {
                    remaining.clear(i);
                    continue;
                }
                if (!system.add(constraints[i].e1, true)) {
                    ++stats.conflicts;
                    return false;
                }
                ++stats.forced;
                remaining.clear(i);
                active = remaining;
                restarted = true;
                break;
            }
            int score = std::popcount(a.remainder | b.remainder);
            if (score < best_score) { best_score = score; branch = i; }
          }
        }
        if (restarted) continue;
        if (remaining.empty()) return true;

        // Disjoint exhaustive partition of (e1=1) OR (e2=1):
        // first e1=1, then e1=0 and e2=1.
        AffineSystem left = system;
        if (left.add(constraints[branch].e1, true) &&
            solve_active(constraints, left, stats, remaining)) return true;
        AffineSystem right = system;
        if (right.add(constraints[branch].e1, false) &&
            right.add(constraints[branch].e2, true) &&
            solve_active(constraints, right, stats, remaining)) return true;
        ++stats.conflicts;
        return false;
    }
}

static bool solve_constraints(const std::vector<Constraint>& constraints,
                              AffineSystem system, SolveStats& stats) {
    ActiveSet active;
    for (unsigned i = 0; i < constraints.size(); ++i)
        active.bits[i / 64] |= 1ull << (i % 64);
    return solve_active(constraints, system, stats, active);
}

static u64 expression(u16 message, unsigned a) {
    u64 e = 0;
    for (int i = 0; i < 12; ++i) if ((message >> i) & 1u)
        for (int j = 0; j < 4; ++j) if ((a >> j) & 1u)
            e |= 1ull << (4 * i + j);
    return e;
}

static bool brute_small(const std::vector<Constraint>& cs, int variables) {
    for (u64 x = 0; x < (1ull << variables); ++x) {
        bool ok = true;
        for (auto c : cs) {
            bool a = std::popcount(x & c.e1) & 1u;
            bool b = std::popcount(x & c.e2) & 1u;
            if (!a && !b) { ok = false; break; }
        }
        if (ok) return true;
    }
    return false;
}

static void self_test() {
    u64 state = 0x9e3779b97f4a7c15ull;
    auto rnd = [&]() {
        state ^= state << 7; state ^= state >> 9; state ^= state << 8; return state;
    };
    for (int trial = 0; trial < 300; ++trial) {
        int n = 1 + int(rnd() % 10);
        int m = int(rnd() % 25);
        std::vector<Constraint> cs;
        u64 mask = (1ull << n) - 1;
        for (int i = 0; i < m; ++i) {
            u64 a = rnd() & mask, b = rnd() & mask;
            if (!a) a = 1ull << (rnd() % n);
            if (!b) b = 1ull << (rnd() % n);
            cs.push_back({a, b});
        }
        SolveStats st;
        bool got = solve_constraints(cs, AffineSystem{}, st);
        bool want = brute_small(cs, n);
        if (got != want) {
            std::cerr << "self-test mismatch on trial " << trial << "\n";
            std::exit(4);
        }
    }
}

static std::string hex4(u16 x) {
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(4) << std::setfill('0') << x;
    return os.str();
}

struct Totals {
    u64 processed = 0, sat = 0, unsat = 0, constraints_sum = 0;
    u64 nodes = 0, conflicts = 0, forced = 0;
    unsigned min_constraints = std::numeric_limits<unsigned>::max();
    unsigned max_constraints = 0;
};

int main(int argc, char** argv) {
    int cls = 2;
    std::string mode = "direct";
    u64 start = 0, end = 2794155;
    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--class" && i + 1 < argc) cls = std::atoi(argv[++i]);
        else if (a == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (a == "--start" && i + 1 < argc) start = std::strtoull(argv[++i], nullptr, 10);
        else if (a == "--end" && i + 1 < argc) end = std::strtoull(argv[++i], nullptr, 10);
        else if (a == "--threads" && i + 1 < argc) threads = std::max(1, std::atoi(argv[++i]));
        else { std::cerr << "usage: --class 1|2|3 --mode direct|orbits --start N --end N --threads N\n"; return 2; }
    }
    if (cls < 1 || cls > 3) { std::cerr << "class must be 1, 2, or 3\n"; return 2; }
    if (mode != "direct" && mode != "orbits") { std::cerr << "unknown mode\n"; return 2; }
    self_test();

    const std::vector<u16> generators = cls == 1
        ? std::vector<u16>{0x8421, 0x3842, 0x6384, 0xc638}
        : cls == 2 ? std::vector<u16>{0x2ba7, 0x8421, 0xe395, 0x16fb}
                   : std::vector<u16>{0x263e, 0x53b7, 0x8421, 0x3842};
    auto basis = dual_basis(generators);
    if (basis.size() != 12 || binary_rank(basis, 16) != 12) {
        std::cerr << "dual basis dimension failure\n"; return 3;
    }
    std::array<u64, 5> rank_distribution{};
    std::vector<RankTwoWord> rank_two;
    for (unsigned x = 0; x < 4096; ++x) {
        u16 m = encode(x, basis);
        int r = rank4(m);
        ++rank_distribution[r];
        if (r == 1) { std::cerr << "dual code contains rank-one word\n"; return 3; }
        if (r == 2) rank_two.push_back(describe_rank_two(x, m));
    }

    auto kernels = enumerate_kernels();
    if (kernels.size() != 2794155) { std::cerr << "kernel enumeration count failure\n"; return 3; }
    std::vector<OrbitRepresentative> orbit_reps;
    std::size_t automorphism_count = 0;
    if (mode == "orbits") {
        auto autos = spread_automorphisms(generators);
        automorphism_count = autos.size();
        auto actions = functional_actions(autos, basis);
        orbit_reps = kernel_orbits(kernels, actions);
        std::cerr << "automorphisms=" << automorphism_count << " kernel_orbits="
                  << orbit_reps.size() << " mass=" << kernels.size() << "\n";
    }
    const u64 item_count = mode == "direct" ? kernels.size() : orbit_reps.size();
    end = std::min<u64>(end, item_count);
    if (start > end) { std::cerr << "invalid range\n"; return 2; }
    const u64 count = end - start;
    std::vector<std::uint16_t> results(count, 0xffffu);
    std::atomic<u64> next{start}, done{0};
    std::atomic<bool> found_sat{false};
    std::mutex totals_mutex, print_mutex;
    Totals global;
    const auto began = std::chrono::steady_clock::now();

    auto worker = [&](unsigned tid) {
        Totals local;
        std::vector<Constraint> constraints;
        constraints.reserve(rank_two.size());
        for (;;) {
            u64 lo = next.fetch_add(64);
            if (lo >= end) break;
            u64 hi = std::min<u64>(lo + 64, end);
            for (u64 index = lo; index < hi; ++index) {
                const u64 kernel_index = mode == "direct" ? index : orbit_reps[index].index;
                const Kernel k = kernels[kernel_index];
                constraints.clear();
                for (const auto& w : rank_two) {
                    if ((std::popcount(u16(w.message & k.h1)) & 1u) ||
                        (std::popcount(u16(w.message & k.h2)) & 1u)) continue;
                    constraints.push_back({expression(w.message, w.annihilator_a),
                                           expression(w.message, w.annihilator_b)});
                }
                SolveStats st;
                bool sat = solve_constraints(constraints, AffineSystem{}, st);
                results[index - start] = static_cast<std::uint16_t>(constraints.size() | (sat ? 0x8000u : 0));
                ++local.processed;
                local.sat += sat;
                local.unsat += !sat;
                local.constraints_sum += constraints.size();
                local.min_constraints = std::min<unsigned>(local.min_constraints, constraints.size());
                local.max_constraints = std::max<unsigned>(local.max_constraints, constraints.size());
                local.nodes += st.nodes; local.conflicts += st.conflicts; local.forced += st.forced;
                if (sat) {
                    found_sat = true;
                    std::lock_guard<std::mutex> lock(print_mutex);
                    std::cerr << "SAT item=" << index << " kernel index=" << kernel_index
                              << " h1=" << hex4(k.h1)
                              << " h2=" << hex4(k.h2) << " constraints=" << constraints.size() << "\n";
                }
            }
            u64 d = done.fetch_add(hi - lo) + (hi - lo);
            if (tid == 0 && d / 100000 != (d - (hi - lo)) / 100000) {
                std::lock_guard<std::mutex> lock(print_mutex);
                double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
                std::cerr << "progress " << d << "/" << count << " elapsed_s=" << sec << "\n";
            }
        }
        std::lock_guard<std::mutex> lock(totals_mutex);
        global.processed += local.processed; global.sat += local.sat; global.unsat += local.unsat;
        global.constraints_sum += local.constraints_sum; global.nodes += local.nodes;
        global.conflicts += local.conflicts; global.forced += local.forced;
        global.min_constraints = std::min(global.min_constraints, local.min_constraints);
        global.max_constraints = std::max(global.max_constraints, local.max_constraints);
    };

    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& t : pool) t.join();
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();

    // Deterministic FNV-1a digest of every covered
    // (item index,kernel index,h1,h2,result,orbit mass) record.
    u64 digest = 14695981039346656037ull;
    auto mix_byte = [&](unsigned char c) { digest ^= c; digest *= 1099511628211ull; };
    for (u64 index = start; index < end; ++index) {
        u64 kernel_index = mode == "direct" ? index : orbit_reps[index].index;
        u64 mass = mode == "direct" ? 1 : orbit_reps[index].mass;
        u64 fields[3] = {index, kernel_index, u64(kernels[kernel_index].h1)
                               | (u64(kernels[kernel_index].h2) << 16)
                               | (u64(results[index - start]) << 32) | (mass << 48)};
        for (u64 f : fields) for (int b = 0; b < 8; ++b) mix_byte((f >> (8 * b)) & 255u);
    }

    u64 coverage_mass = 0;
    for (u64 index = start; index < end; ++index)
        coverage_mass += mode == "direct" ? 1 : orbit_reps[index].mass;
    std::cout << "{\"checker\":\"affine-equation-dpll-48var-v3\",\"mode\":\""
              << mode << "\",\"class\":" << cls
              << ",\"generators\":[";
    for (std::size_t i = 0; i < generators.size(); ++i)
        std::cout << (i ? "," : "") << "\"" << hex4(generators[i]) << "\"";
    std::cout << "],\"dual_basis\":[";
    for (std::size_t i = 0; i < basis.size(); ++i)
        std::cout << (i ? "," : "") << "\"" << hex4(basis[i]) << "\"";
    std::cout << "],\"dual_rank_distribution\":[";
    for (int i = 0; i < 5; ++i) std::cout << (i ? "," : "") << rank_distribution[i];
    std::cout << "],\"rank_two_words\":" << rank_two.size()
              << ",\"total_kernels\":" << kernels.size()
              << ",\"automorphism_count\":" << automorphism_count
              << ",\"total_orbits\":" << (mode == "orbits" ? orbit_reps.size() : 0)
              << ",\"range_start\":" << start
              << ",\"range_end\":" << end << ",\"processed\":" << global.processed
              << ",\"coverage_mass\":" << coverage_mass
              << ",\"unsat\":" << global.unsat << ",\"sat\":" << global.sat
              << ",\"constraint_count_min\":" << global.min_constraints
              << ",\"constraint_count_max\":" << global.max_constraints
              << ",\"constraint_count_sum\":" << global.constraints_sum
              << ",\"search_nodes\":" << global.nodes << ",\"conflicts\":" << global.conflicts
              << ",\"forced_equations\":" << global.forced << ",\"threads\":" << threads
              << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(6) << elapsed
              << ",\"record_digest_fnv1a64\":\"0x" << std::hex << std::setw(16)
              << std::setfill('0') << digest << "\"}\n";
    return found_sat ? 1 : 0;
}
