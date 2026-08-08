#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using u16 = std::uint16_t;
using u64 = std::uint64_t;
static constexpr u16 IDENTITY = 0x8421;

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

static u16 multiply4(u16 a, u16 b) {
    u16 z = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned ar = (a >> (4 * i)) & 15u, zr = 0;
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

static u16 transpose4(u16 a) {
    u16 z = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            z |= ((a >> (4 * i + j)) & 1u) << (4 * j + i);
    return z;
}

static std::vector<u16> span(const std::vector<u16>& basis) {
    std::vector<u16> v(1u << basis.size());
    for (unsigned x = 0; x < v.size(); ++x) {
        u16 a = 0;
        for (unsigned i = 0; i < basis.size(); ++i) if ((x >> i) & 1u) a ^= basis[i];
        v[x] = a;
    }
    return v;
}

// Linear section of M_4/<I>: clear bit zero by adding I when necessary.
static u16 quotient(u16 a) { return (a & 1u) ? u16(a ^ IDENTITY) : a; }

struct Spread { u16 a, b, c; };
static u64 key(Spread s) { return u64(s.a) | (u64(s.b) << 16) | (u64(s.c) << 32); }

static Spread canonical_spread(const std::vector<u16>& generators) {
    auto words = span(generators);
    std::vector<u16> q;
    for (u16 a : words) {
        u16 x = quotient(a);
        if (x) q.push_back(x);
    }
    std::sort(q.begin(), q.end());
    q.erase(std::unique(q.begin(), q.end()), q.end());
    if (q.size() != 7) { std::cerr << "bad quotient dimension\n"; std::exit(3); }
    u16 a = q[0], b = q[1], d = a ^ b, c = 0;
    for (u16 x : q) if (x != a && x != b && x != d) { c = x; break; }
    if (!c) { std::cerr << "canonical basis failure\n"; std::exit(3); }
    return {a, b, c};
}

static std::string hex4(u16 x) {
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(4) << std::setfill('0') << x;
    return os.str();
}

int main() {
    std::array<bool, 65536> good{};
    std::array<int, 65536> good_index{};
    good_index.fill(-1);
    std::vector<u16> vertices;
    for (unsigned x = 0; x < 65536; x += 2) {
        if (rank4(x) == 4 && rank4(u16(x ^ IDENTITY)) == 4) {
            good[x] = true;
            good_index[x] = vertices.size();
            vertices.push_back(static_cast<u16>(x));
        }
    }
    const std::size_t words = (vertices.size() + 63) / 64;
    std::vector<u64> adjacency(vertices.size() * words);
    for (std::size_t i = 0; i < vertices.size(); ++i)
        for (std::size_t j = 0; j < vertices.size(); ++j)
            if (good[vertices[i] ^ vertices[j]])
                adjacency[i * words + j / 64] |= 1ull << (j % 64);

    std::vector<Spread> spreads;
    for (std::size_t ia = 0; ia < vertices.size(); ++ia) {
        u16 a = vertices[ia];
        for (std::size_t ib = ia + 1; ib < vertices.size(); ++ib) {
            u16 b = vertices[ib], d = a ^ b;
            if (!good[d] || d <= b) continue; // a,b are the two smallest in their plane.
            int id = good_index[d];
            for (std::size_t w = 0; w < words; ++w) {
                u64 candidates = adjacency[ia * words + w] & adjacency[ib * words + w]
                               & adjacency[id * words + w];
                while (candidates) {
                    int bit = std::countr_zero(candidates);
                    candidates &= candidates - 1;
                    std::size_t ic = 64 * w + bit;
                    if (ic >= vertices.size()) continue;
                    u16 c = vertices[ic];
                    if (c <= b) continue;
                    if (c < (c ^ a) && c < (c ^ b) && c < (c ^ d))
                        spreads.push_back({a, b, c});
                }
            }
        }
    }
    std::sort(spreads.begin(), spreads.end(), [](Spread x, Spread y) { return key(x) < key(y); });
    spreads.erase(std::unique(spreads.begin(), spreads.end(),
                             [](Spread x, Spread y) { return key(x) == key(y); }), spreads.end());

    std::unordered_map<u64, std::size_t> index;
    for (std::size_t i = 0; i < spreads.size(); ++i) index.emplace(key(spreads[i]), i);
    std::vector<u16> gl;
    for (unsigned a = 0; a < 65536; ++a) if (rank4(a) == 4) gl.push_back(a);

    std::vector<bool> visited(spreads.size());
    std::vector<std::vector<std::size_t>> orbits;
    for (std::size_t root = 0; root < spreads.size(); ++root) if (!visited[root]) {
        Spread s = spreads[root];
        std::vector<u16> basis{IDENTITY, s.a, s.b, s.c};
        auto elements = span(basis);
        std::vector<std::size_t> orbit;
        for (u16 a : elements) if (a) {
            u16 ai = *inverse4(a);
            for (u16 q : gl) {
                u16 p = multiply4(*inverse4(q), ai);
                std::vector<u16> image;
                for (u16 g : basis) image.push_back(multiply4(multiply4(p, g), q));
                auto it = index.find(key(canonical_spread(image)));
                if (it == index.end()) { std::cerr << "orbit image missing\n"; return 3; }
                orbit.push_back(it->second);
            }
        }
        std::sort(orbit.begin(), orbit.end());
        orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
        for (auto j : orbit) {
            if (visited[j]) { std::cerr << "orbit overlap\n"; return 3; }
            visited[j] = true;
        }
        orbits.push_back(std::move(orbit));
    }

    const std::array<std::vector<u16>, 3> named{{
        {0x8421, 0x3842, 0x6384, 0xc638},
        {0x2ba7, 0x8421, 0xe395, 0x16fb},
        {0x263e, 0x53b7, 0x8421, 0x3842}
    }};
    std::cout << "{\"good_quotient_vectors\":" << vertices.size()
              << ",\"normalized_spreads\":" << spreads.size()
              << ",\"lr_orbits\":" << orbits.size() << ",\"orbits\":[";
    for (std::size_t o = 0; o < orbits.size(); ++o) {
        std::size_t r = orbits[o][0];
        auto ts = canonical_spread({IDENTITY, transpose4(spreads[r].a),
                                    transpose4(spreads[r].b), transpose4(spreads[r].c)});
        std::size_t ti = index.at(key(ts)), to = 0;
        while (to < orbits.size() && !std::binary_search(orbits[to].begin(), orbits[to].end(), ti)) ++to;
        std::cout << (o ? "," : "") << "{\"root_index\":" << r
                  << ",\"mass\":" << orbits[o].size() << ",\"transpose_orbit\":" << to << "}";
    }
    std::cout << "],\"named_classes\":[";
    for (std::size_t n = 0; n < named.size(); ++n) {
        std::size_t i = index.at(key(canonical_spread(named[n]))), o = 0;
        while (o < orbits.size() && !std::binary_search(orbits[o].begin(), orbits[o].end(), i)) ++o;
        std::cout << (n ? "," : "") << "{\"name\":" << (n + 1)
                  << ",\"normalized_index\":" << i << ",\"orbit\":" << o << "}";
    }
    std::cout << "]}\n";
    return 0;
}
