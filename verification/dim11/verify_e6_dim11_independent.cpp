#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::size_t value_start(const std::string& text, const std::string& key) {
  const std::string quoted = "\"" + key + "\"";
  std::size_t pos = text.find(quoted);
  if (pos == std::string::npos) throw std::runtime_error("missing key: " + key);
  pos = text.find(':', pos + quoted.size());
  if (pos == std::string::npos) throw std::runtime_error("missing colon: " + key);
  return pos + 1;
}

void skip_space(const std::string& text, std::size_t& pos) {
  while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
}

unsigned parse_unsigned(const std::string& text, std::size_t& pos) {
  skip_space(text, pos);
  if (pos == text.size() || !std::isdigit(static_cast<unsigned char>(text[pos])))
    throw std::runtime_error("expected unsigned integer");
  unsigned result = 0;
  while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) {
    result = 10 * result + unsigned(text[pos] - '0');
    ++pos;
  }
  return result;
}

std::vector<unsigned> parse_vector(const std::string& text, const std::string& key) {
  std::size_t pos = value_start(text, key);
  skip_space(text, pos);
  if (pos == text.size() || text[pos++] != '[') throw std::runtime_error("expected vector");
  std::vector<unsigned> result;
  for (;;) {
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == ']') { ++pos; return result; }
    result.push_back(parse_unsigned(text, pos));
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == ',') { ++pos; continue; }
    if (pos < text.size() && text[pos] == ']') { ++pos; return result; }
    throw std::runtime_error("malformed vector");
  }
}

std::vector<std::vector<unsigned>> parse_matrix(const std::string& text,
                                                const std::string& key) {
  std::size_t pos = value_start(text, key);
  skip_space(text, pos);
  if (pos == text.size() || text[pos++] != '[') throw std::runtime_error("expected matrix");
  std::vector<std::vector<unsigned>> result;
  for (;;) {
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == ']') { ++pos; return result; }
    if (pos == text.size() || text[pos++] != '[') throw std::runtime_error("expected row");
    std::vector<unsigned> row;
    for (;;) {
      skip_space(text, pos);
      if (pos < text.size() && text[pos] == ']') { ++pos; break; }
      row.push_back(parse_unsigned(text, pos));
      skip_space(text, pos);
      if (pos < text.size() && text[pos] == ',') { ++pos; continue; }
      if (pos < text.size() && text[pos] == ']') { ++pos; break; }
      throw std::runtime_error("malformed row");
    }
    result.push_back(std::move(row));
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == ',') { ++pos; continue; }
    if (pos < text.size() && text[pos] == ']') { ++pos; return result; }
    throw std::runtime_error("malformed matrix");
  }
}

unsigned parse_scalar(const std::string& text, const std::string& key) {
  std::size_t pos = value_start(text, key);
  return parse_unsigned(text, pos);
}

std::map<unsigned, unsigned> parse_distribution(const std::string& text) {
  std::size_t pos = value_start(text, "rank_distribution");
  skip_space(text, pos);
  if (pos == text.size() || text[pos++] != '{') throw std::runtime_error("expected object");
  std::map<unsigned, unsigned> result;
  for (;;) {
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == '}') { ++pos; return result; }
    if (pos == text.size() || text[pos++] != '"') throw std::runtime_error("expected rank key");
    unsigned rank = parse_unsigned(text, pos);
    if (pos == text.size() || text[pos++] != '"') throw std::runtime_error("bad rank key");
    skip_space(text, pos);
    if (pos == text.size() || text[pos++] != ':') throw std::runtime_error("missing rank colon");
    unsigned count = parse_unsigned(text, pos);
    if (!result.emplace(rank, count).second) throw std::runtime_error("duplicate rank");
    skip_space(text, pos);
    if (pos < text.size() && text[pos] == ',') { ++pos; continue; }
    if (pos < text.size() && text[pos] == '}') { ++pos; return result; }
    throw std::runtime_error("malformed distribution");
  }
}

unsigned binary_rank(std::vector<uint32_t> rows, int bit_count) {
  unsigned rank = 0;
  for (int bit = bit_count - 1; bit >= 0; --bit) {
    std::size_t pivot = rank;
    while (pivot < rows.size() && ((rows[pivot] >> bit) & 1U) == 0) ++pivot;
    if (pivot == rows.size()) continue;
    std::swap(rows[rank], rows[pivot]);
    for (std::size_t i = 0; i < rows.size(); ++i)
      if (i != rank && ((rows[i] >> bit) & 1U)) rows[i] ^= rows[rank];
    ++rank;
  }
  return rank;
}

unsigned matrix_rank(const std::vector<unsigned>& columns, unsigned row_count) {
  std::vector<uint32_t> rows(row_count, 0);
  for (std::size_t column = 0; column < columns.size(); ++column)
    for (unsigned row = 0; row < row_count; ++row)
      if ((columns[column] >> row) & 1U) rows[row] |= uint32_t(1) << column;
  return binary_rank(std::move(rows), int(columns.size()));
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::runtime_error("usage: verifier certificate.json");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open certificate");
    const std::string text((std::istreambuf_iterator<char>(input)), {});

    const auto basis = parse_matrix(text, "basis_column_masks");
    const auto heights = parse_vector(text, "support_column_heights");
    const unsigned claimed_dimension = parse_scalar(text, "dimension");
    const unsigned claimed_minimum = parse_scalar(text, "minimum_nonzero_rank");
    const auto claimed_distribution = parse_distribution(text);

    if (basis.size() != claimed_dimension) throw std::runtime_error("dimension mismatch");
    if (heights.empty()) throw std::runtime_error("empty support");
    const unsigned row_count = *std::max_element(heights.begin(), heights.end());

    std::vector<uint32_t> packed_basis;
    for (const auto& generator : basis) {
      if (generator.size() != heights.size()) throw std::runtime_error("generator width mismatch");
      uint32_t packed = 0;
      for (std::size_t column = 0; column < generator.size(); ++column) {
        if (generator[column] >= (uint32_t(1) << heights[column]))
          throw std::runtime_error("generator violates Ferrers support");
        packed |= uint32_t(generator[column]) << (row_count * column);
      }
      packed_basis.push_back(packed);
    }
    const unsigned generator_rank = binary_rank(packed_basis, int(row_count * heights.size()));
    if (generator_rank != claimed_dimension) throw std::runtime_error("dependent generators");

    std::set<uint32_t> distinct_words;
    std::map<unsigned, unsigned> distribution;
    unsigned observed_minimum = row_count + 1;
    const uint32_t word_count = uint32_t(1) << claimed_dimension;
    for (uint32_t coefficients = 0; coefficients < word_count; ++coefficients) {
      std::vector<unsigned> columns(heights.size(), 0);
      uint32_t packed = 0;
      for (unsigned i = 0; i < basis.size(); ++i) if ((coefficients >> i) & 1U) {
        packed ^= packed_basis[i];
        for (std::size_t column = 0; column < columns.size(); ++column)
          columns[column] ^= basis[i][column];
      }
      distinct_words.insert(packed);
      const unsigned rank = matrix_rank(columns, row_count);
      ++distribution[rank];
      if (coefficients != 0) observed_minimum = std::min(observed_minimum, rank);
    }

    if (distinct_words.size() != word_count) throw std::runtime_error("word collision");
    if (observed_minimum != claimed_minimum) throw std::runtime_error("minimum-rank mismatch");
    if (distribution != claimed_distribution) throw std::runtime_error("distribution mismatch");

    std::cout << "{\"status\":\"VERIFIED\",\"support_ok\":true,"
              << "\"generator_rank\":" << generator_rank
              << ",\"distinct_codewords\":" << distinct_words.size()
              << ",\"minimum_nonzero_rank\":" << observed_minimum
              << ",\"rank_distribution\":{";
    bool first = true;
    for (const auto& [rank, count] : distribution) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '\"' << rank << "\":" << count;
    }
    std::cout << "}}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "verification failed: " << error.what() << '\n';
    return 1;
  }
}
