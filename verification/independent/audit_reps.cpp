#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

using u16 = std::uint16_t;

static int rank4(u16 a) {
    std::array<unsigned, 4> r{};
    for (int i = 0; i < 4; ++i) r[i] = (a >> (4 * i)) & 15u;
    int rank = 0;
    for (int c = 0; c < 4; ++c) {
        int p = rank;
        while (p < 4 && ((r[p] >> c) & 1u) == 0) ++p;
        if (p == 4) continue;
        std::swap(r[p], r[rank]);
        for (int i = 0; i < 4; ++i)
            if (i != rank && ((r[i] >> c) & 1u)) r[i] ^= r[rank];
        ++rank;
    }
    return rank;
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
    for (int i = 0; i < 4; ++i)
        r[i] = ((a >> (4 * i)) & 15u) | (1u << (4 + i));
    for (int c = 0; c < 4; ++c) {
        int p = c;
        while (p < 4 && ((r[p] >> c) & 1u) == 0) ++p;
        if (p == 4) return std::nullopt;
        std::swap(r[p], r[c]);
        for (int i = 0; i < 4; ++i)
            if (i != c && ((r[i] >> c) & 1u)) r[i] ^= r[c];
    }
    u16 z = 0;
    for (int i = 0; i < 4; ++i) z |= ((r[i] >> 4) & 15u) << (4 * i);
    return z;
}

static std::vector<u16> span(const std::vector<u16>& basis) {
    std::vector<u16> v(1u << basis.size());
    for (std::size_t x = 0; x < v.size(); ++x) {
        u16 a = 0;
        for (std::size_t i = 0; i < basis.size(); ++i)
            if ((x >> i) & 1u) a ^= basis[i];
        v[x] = a;
    }
    std::sort(v.begin(), v.end());
    return v;
}

static int binary_rank(std::vector<u16> rows, int width) {
    int rank = 0;
    for (int c = 0; c < width; ++c) {
        int p = rank;
        while (p < static_cast<int>(rows.size()) && ((rows[p] >> c) & 1u) == 0) ++p;
        if (p == static_cast<int>(rows.size())) continue;
        std::swap(rows[p], rows[rank]);
        for (int i = 0; i < static_cast<int>(rows.size()); ++i)
            if (i != rank && ((rows[i] >> c) & 1u)) rows[i] ^= rows[rank];
        ++rank;
    }
    return rank;
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
    std::array<bool, 16> is_pivot{};
    for (int p : pivots) is_pivot[p] = true;
    std::vector<u16> basis;
    for (int f = 0; f < 16; ++f) {
        if (is_pivot[f]) continue;
        u16 x = u16(1u << f);
        for (int i = 0; i < r; ++i)
            if ((equations[i] >> f) & 1u) x |= u16(1u << pivots[i]);
        basis.push_back(x);
    }
    return basis;
}

static std::string hex4(u16 x) {
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(4) << std::setfill('0') << x;
    return os.str();
}

static bool same_set(std::vector<u16> a, std::vector<u16> b) {
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
}

struct LRMap { u16 p, q; };

static std::optional<LRMap> find_lr_map(const std::vector<u16>& from,
                                        const std::vector<u16>& to) {
    const auto from_set = span(from);
    const auto to_set = span(to);
    const u16 identity = 0x8421;
    for (u16 q = 0; ; ++q) {
        auto qi = inverse4(q);
        if (qi) {
            for (u16 b : to_set) {
                if (b == 0) continue;
                // Force P I Q = b, so P = b Q^{-1}.
                u16 p = multiply4(b, *qi);
                std::vector<u16> image;
                image.reserve(16);
                for (u16 a : from_set) image.push_back(multiply4(multiply4(p, a), q));
                if (same_set(image, to_set)) return LRMap{p, q};
            }
        }
        if (q == 0xffff) break;
    }
    (void)identity;
    return std::nullopt;
}

static void audit(const std::string& name, const std::vector<u16>& generators) {
    auto s = span(generators);
    auto dual = dual_basis(generators);
    auto u = span(dual);
    std::array<std::uint64_t, 5> sranks{}, uranks{};
    bool orthogonal = true;
    for (u16 a : s) ++sranks[rank4(a)];
    for (u16 a : u) {
        ++uranks[rank4(a)];
        for (u16 g : generators)
            if (std::popcount(u16(a & g)) & 1u) orthogonal = false;
    }
    std::vector<u16> transposed;
    for (u16 g : generators) transposed.push_back(transpose4(g));
    auto self_t = find_lr_map(transposed, generators);

    std::cout << "{\"name\":\"" << name << "\",\"generator_rank\":"
              << binary_rank(generators, 16) << ",\"spread_rank_distribution\":[";
    for (int i = 0; i < 5; ++i) std::cout << (i ? "," : "") << sranks[i];
    std::cout << "],\"dual_basis\":[";
    for (std::size_t i = 0; i < dual.size(); ++i)
        std::cout << (i ? "," : "") << "\"" << hex4(dual[i]) << "\"";
    std::cout << "],\"dual_dimension\":" << binary_rank(dual, 16)
              << ",\"dual_rank_distribution\":[";
    for (int i = 0; i < 5; ++i) std::cout << (i ? "," : "") << uranks[i];
    std::cout << "],\"orthogonality_ok\":" << (orthogonal ? "true" : "false")
              << ",\"transpose_lr_equivalent\":" << (self_t ? "true" : "false");
    if (self_t)
        std::cout << ",\"transpose_map\":{\"P\":\"" << hex4(self_t->p)
                  << "\",\"Q\":\"" << hex4(self_t->q) << "\"}";
    std::cout << "}\n";
}

int main() {
    const std::vector<u16> s1{0x8421, 0x3842, 0x6384, 0xc638};
    const std::vector<u16> s2{0x2ba7, 0x8421, 0xe395, 0x16fb};
    const std::vector<u16> s3{0x263e, 0x53b7, 0x8421, 0x3842};
    audit("class1_field", s1);
    audit("class2", s2);
    audit("class3", s3);
    auto eq12 = find_lr_map(s1, s2);
    auto eq13 = find_lr_map(s1, s3);
    auto eq23 = find_lr_map(s2, s3);
    std::cout << "{\"class1_lr_equivalent_class2\":" << (eq12 ? "true" : "false")
              << ",\"class1_lr_equivalent_class3\":" << (eq13 ? "true" : "false") << "}\n";
    std::cout << "{\"class2_lr_equivalent_class3\":" << (eq23 ? "true" : "false") << "}\n";
    return 0;
}
