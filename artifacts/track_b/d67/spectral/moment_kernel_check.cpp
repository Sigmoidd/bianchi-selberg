// Exact sufficient check for the constant local-CR energy kernel.
// Run this only on independently verified positive dyadic moment rows.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>
using U = uint32_t;
const U NIL = UINT32_MAX;
template <class T> void read(std::ifstream &f, T &x) {
  f.read(reinterpret_cast<char *>(&x), sizeof(x));
  if (!f) throw std::runtime_error("truncated moment rows");
}
struct DSU {
  std::vector<U> p, size;
  U count;
  DSU(U n) : p(n), size(n, 1), count(n) {
    for (U i = 0; i < n; i++) p[i] = i;
  }
  U root(U i) {
    while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
    return i;
  }
  bool join(U i, U j) {
    i = root(i); j = root(j);
    if (i == j) return false;
    if (size[i] < size[j]) std::swap(i, j);
    p[j] = i; size[i] += size[j]; count--; return true;
  }
};
struct Row { U tet; std::vector<U> masters; };
int main(int argc, char **argv) {
  try {
    if (argc != 2) throw std::runtime_error("usage: kernel_check ROWS");
    std::ifstream f(argv[1], std::ios::binary);
    uint64_t nr; read(f, nr);
    if (nr % 4 || nr > 40000000) throw std::runtime_error("invalid row count");
    DSU d(nr / 4); std::vector<U> owners; std::vector<Row> rows;
    for (uint64_t i = 0; i < nr; i++) {
      U n; read(f, n); Row row{U(i / 4), {}};
      if (!n || n > 1000000) throw std::runtime_error("invalid support");
      row.masters.reserve(n);
      for (U j = 0; j < n; j++) {
        U c, e; uint64_t num; read(f, c); read(f, num); read(f, e);
        if (!num || e > 48 || c >= nr) throw std::runtime_error("invalid positive row");
        if (c >= owners.size()) owners.resize(c + 1, NIL);
        row.masters.push_back(c);
        if (n == 1) {
          if (num != 1 || e != 0) throw std::runtime_error("identity weight mismatch");
          if (owners[c] == NIL) owners[c] = i / 4;
          else d.join(i / 4, owners[c]);
        }
      }
      if (n > 1) rows.push_back(std::move(row));
    }
    if (f.peek() != EOF) throw std::runtime_error("trailing rows");
    for (U owner : owners)
      if (owner == NIL) throw std::runtime_error("master identity witness missing");
    U initial = d.count, passes = 0; bool change = true;
    while (change) {
      change = false; passes++;
      for (auto &row : rows) {
        U head = d.root(owners[row.masters[0]]); bool same = true;
        for (U c : row.masters)
          if (d.root(owners[c]) != head) { same = false; break; }
        if (same) change = d.join(row.tet, head) || change;
      }
    }
    std::cout << "{\n  \"schema\": \"d67-exact-moment-energy-kernel-check/v1\",\n"
      << "  \"tetrahedra\": " << nr/4 << ",\n"
      << "  \"identity_components\": " << initial << ",\n"
      << "  \"components_after_positive_row_merging\": " << d.count << ",\n"
      << "  \"passes\": " << passes << ",\n"
      << "  \"constant_energy_kernel_proved\": " << (d.count == 1 ? "true" : "false") << ",\n"
      << "  \"spectral_exclusion_certified\": false\n}\n";
    return 0;
  } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
