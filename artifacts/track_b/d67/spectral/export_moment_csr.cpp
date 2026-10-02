// Convert independently verified exact dyadic rows to diagnostic binary CSR.
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
template <class T> void get(std::ifstream &f, T &x) {
  f.read(reinterpret_cast<char *>(&x), sizeof(x));
  if (!f)
    throw std::runtime_error("truncated rows");
}
template <class T> void put(std::ofstream &f, const T &x) {
  f.write(reinterpret_cast<const char *>(&x), sizeof(x));
}
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("usage: export ROWS CSR");
    std::ifstream in(argv[1], std::ios::binary);
    uint64_t rows;
    get(in, rows);
    std::vector<uint64_t> ptr(rows + 1);
    std::vector<uint32_t> col;
    std::vector<double> val;
    uint64_t ndof = 0;
    for (uint64_t r = 0; r < rows; r++) {
      uint32_t n;
      get(in, n);
      for (uint32_t j = 0; j < n; j++) {
        uint32_t c, e;
        uint64_t num;
        get(in, c);
        get(in, num);
        get(in, e);
        if (e > 48 || num > (uint64_t(1) << 53))
          throw std::runtime_error("non-exact diagnostic double weight");
        col.push_back(c);
        val.push_back(std::ldexp(double(num), -int(e)));
        ndof = std::max(ndof, uint64_t(c) + 1);
      }
      ptr[r + 1] = col.size();
    }
    if (in.peek() != EOF)
      throw std::runtime_error("trailing rows");
    std::ofstream out(argv[2], std::ios::binary);
    put(out, rows);
    put(out, ndof);
    uint64_t nz = col.size();
    put(out, nz);
    out.write(reinterpret_cast<const char *>(ptr.data()), ptr.size() * 8);
    out.write(reinterpret_cast<const char *>(col.data()), col.size() * 4);
    out.write(reinterpret_cast<const char *>(val.data()), val.size() * 8);
    if (!out)
      throw std::runtime_error("CSR write failed");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
