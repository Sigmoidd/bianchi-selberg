// Independent exact topology/moment replay: geometric barycentrics, leaf-path
// reconstruction, physical incidence and exact dyadic row checks.
// Build: g++ -std=c++17 -O2 verify_moment_map.cpp -o /tmp/verify-d67-moments
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using I = __int128_t;
using U = uint32_t;
const U NIL = UINT32_MAX;
using Tri = array<U, 3>;
using Tet = array<U, 4>;
void check(bool x, const char *s) {
  if (!x)
    throw runtime_error(s);
}
I decimal(const string &s) {
  I x = 0;
  for (char c : s)
    if (c != '-')
      x = 10 * x + c - '0';
  return s[0] == '-' ? -x : x;
}
template <class T> void read(ifstream &f, T &x) {
  f.read(reinterpret_cast<char *>(&x), sizeof(x));
  check(bool(f), "truncated binary");
}
struct Point {
  array<I, 3> x;
  bool operator==(const Point &p) const { return x == p.x; }
};
uint64_t hashword(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  return x ^ (x >> 33);
}
struct PH {
  size_t operator()(const Point &p) const {
    uint64_t h = 0;
    for (I x : p.x) {
      h = hashword(h ^ uint64_t(x));
      h = hashword(h ^ uint64_t(x >> 64));
    }
    return h;
  }
};
struct FH {
  size_t operator()(const Tri &t) const {
    return hashword(t[0]) ^ hashword(uint64_t(t[1]) + 7777777) ^
           hashword(uint64_t(t[2]) + 999999999);
  }
};
struct Face {
  Tri v;
  U parent;
  array<U, 4> child;
  U cell;
};
struct Key {
  U root, depth;
  array<U, 9> b;
  bool operator==(const Key &k) const {
    return root == k.root && depth == k.depth && b == k.b;
  }
  bool operator<(const Key &k) const {
    return tie(root, depth, b) < tie(k.root, k.depth, k.b);
  }
};
struct Cell {
  Key key;
  array<U, 4> child;
  U dof;
};
struct Seed {
  U rep;
  array<int, 3> frame;
  vector<array<int, 3>> group;
};
using Weight = pair<uint64_t, U>;
void sumweight(Weight &x, Weight y) {
  U e = max(x.second, y.second);
  check(e <= 48, "dyadic exponent overflow");
  uint64_t a = x.first << (e - x.second), b = y.first << (e - y.second);
  check(a <= UINT64_MAX - b, "dyadic sum overflow");
  x = {a + b, e};
  while (x.second && x.first % 2 == 0) {
    x.first /= 2;
    x.second--;
  }
}
vector<Cell> cells;
void expectedrow(U c, U exponent, map<U, Weight> &row, set<U> &active) {
  check(active.insert(c).second, "cyclic canonical hierarchy");
  if (cells[c].child[0] == NIL)
    sumweight(row[cells[c].dof], {1, exponent});
  else
    for (U j : cells[c].child)
      expectedrow(j, exponent + 2, row, active);
  active.erase(c);
}
int main(int argc, char **argv) {
  try {
    check(argc == 3, "usage: verifier INPUT OUTPUT_PREFIX");
    ifstream src(argv[1]);
    size_t nn, nr, ns, nl;
    string zs, rs;
    src >> nn >> nr >> ns >> nl >> zs >> rs;
    vector<Point> initial(nn);
    for (auto &p : initial) {
      int64_t a, b;
      string r;
      src >> a >> b >> r;
      p.x = {a, b, decimal(r)};
    }
    vector<Tet> roots(nr);
    for (auto &t : roots)
      for (auto &n : t)
        src >> n;
    struct RawSeed {
      Tri source, rep;
      array<int, 3> frame;
      vector<array<int, 3>> group;
    };
    vector<RawSeed> raw(ns);
    for (auto &s : raw) {
      for (auto &n : s.source)
        src >> n;
      for (auto &n : s.rep)
        src >> n;
      for (auto &j : s.frame)
        src >> j;
      size_t ng;
      src >> ng;
      s.group.resize(ng);
      for (auto &g : s.group)
        for (auto &j : g)
          src >> j;
    }
    vector<pair<U, string>> plan(nl);
    for (auto &p : plan) {
      src >> p.first >> p.second;
      if (p.second == "*")
        p.second = "";
    }
    check(bool(src), "truncated source input");
    string prefix = argv[2];
    ifstream topo(prefix + ".topology.bin", ios::binary);
    uint64_t magic, nnode, nface, ncell, nleaf, ndof;
    read(topo, magic);
    check(magic == 0x4436374d4f4d3031ULL, "wrong topology schema");
    read(topo, nnode);
    read(topo, nface);
    read(topo, ncell);
    read(topo, nleaf);
    read(topo, ndof);
    check(nleaf == nl && nnode >= nn && nnode < 10000000 && nface < 10000000 &&
              ncell <= nface && ndof <= ncell,
          "invalid counts");
    vector<Point> points(nnode);
    unordered_map<Point, U, PH> pointindex;
    pointindex.reserve(nnode);
    for (U i = 0; i < nnode; i++) {
      int64_t a, b;
      uint64_t lo, hi;
      array<U, 2> parents;
      read(topo, a);
      read(topo, b);
      read(topo, lo);
      read(topo, hi);
      read(topo, parents);
      points[i].x = {a, b, I((__uint128_t(hi) << 64) | lo)};
      if (i < nn) {
        check(points[i] == initial[i] && parents == array<U, 2>{NIL, NIL},
              "initial point mismatch");
      } else {
        check(parents[0] < i && parents[1] < i && parents[0] != parents[1],
              "invalid midpoint parents");
        for (int k = 0; k < 3; k++)
          check(2 * points[i].x[k] ==
                    points[parents[0]].x[k] + points[parents[1]].x[k],
                "midpoint relation fails");
      }
      check(pointindex.emplace(points[i], i).second, "duplicate point");
    }
    vector<Face> faces(nface);
    unordered_map<Tri, U, FH> faceindex;
    faceindex.reserve(nface);
    for (U i = 0; i < nface; i++) {
      auto &f = faces[i];
      read(topo, f.v);
      read(topo, f.parent);
      read(topo, f.child);
      read(topo, f.cell);
      check(f.v[0] < f.v[1] && f.v[1] < f.v[2] && f.v[2] < nnode &&
                f.cell < ncell,
            "invalid face");
      check(faceindex.emplace(f.v, i).second, "duplicate physical face");
      check(f.parent == NIL || f.parent < nface, "face parent out of range");
      for (U c : f.child)
        check(c == NIL || c < nface, "child out of range");
      check(all_of(f.child.begin(), f.child.end(),
                   [&](U c) { return (c == NIL) == (f.child[0] == NIL); }),
            "partial subdivision");
    }
    auto midpoint = [&](U i, U j) {
      Point p;
      for (int k = 0; k < 3; k++) {
        I sum = points[i].x[k] + points[j].x[k];
        check(sum % 2 == 0, "nonintegral midpoint");
        p.x[k] = sum / 2;
      }
      auto it = pointindex.find(p);
      check(it != pointindex.end(), "missing midpoint");
      return it->second;
    };
    auto lookupface = [&](Tri v) {
      sort(v.begin(), v.end());
      auto it = faceindex.find(v);
      check(it != faceindex.end(), "missing triangle");
      return it->second;
    };
    for (U i = 0; i < nface; i++)
      if (faces[i].child[0] != NIL) {
        auto v = faces[i].v;
        U a = v[0], b = v[1], c = v[2], ab = midpoint(a, b),
          ac = midpoint(a, c), bc = midpoint(b, c);
        array<Tri, 4> q = {Tri{a, ab, ac}, Tri{b, ab, bc}, Tri{c, ac, bc},
                           Tri{ab, ac, bc}};
        set<U> found;
        for (auto t : q) {
          U child = lookupface(t);
          check(faces[child].parent == i, "quarter-face parent mismatch");
          found.insert(child);
        }
        check(found == set<U>(faces[i].child.begin(), faces[i].child.end()),
              "face children fail exact quarter partition");
      }
    for (U i = 0; i < nface; i++)
      if (faces[i].parent != NIL) {
        auto ch = faces[faces[i].parent].child;
        check(find(ch.begin(), ch.end(), i) != ch.end(), "parent omits child");
      }
    cells.resize(ncell);
    set<Key> unique;
    vector<uint8_t> dofs(ndof);
    for (auto &c : cells) {
      read(topo, c.key.root);
      read(topo, c.key.depth);
      read(topo, c.key.b);
      read(topo, c.child);
      read(topo, c.dof);
      check(c.key.root < nface && faces[c.key.root].parent == NIL &&
                c.key.depth <= 16,
            "invalid canonical root");
      check(unique.insert(c.key).second, "duplicate canonical cell");
      for (U j : c.child)
        check(j == NIL || j < ncell, "canonical child out of range");
      if (c.child[0] == NIL) {
        check(c.dof < ndof && !dofs[c.dof], "duplicate or invalid master dof");
        dofs[c.dof] = 1;
      } else {
        check(c.dof == NIL, "nonterminal master dof");
        for (U j : c.child)
          check(cells.size() > j, "bad canonical child");
      }
    }
    check(all_of(dofs.begin(), dofs.end(), [](auto b) { return b == 1; }),
          "missing dof");
    for (auto &c : cells) {
      check(all_of(c.child.begin(), c.child.end(),
                   [&](U j) { return (j == NIL) == (c.child[0] == NIL); }),
            "partial canonical subdivision");
      if (c.child[0] != NIL)
        for (U j : c.child)
          check(cells[j].key.root == c.key.root &&
                    cells[j].key.depth == c.key.depth + 1,
                "invalid canonical child depth");
    }
    unordered_map<U, Seed> seeds;
    for (auto &s : raw) {
      U f = lookupface(s.source), rep = lookupface(s.rep);
      check(faces[f].parent == NIL && faces[rep].parent == NIL, "nonroot seed");
      check(seeds.emplace(f, Seed{rep, s.frame, s.group}).second,
            "duplicate seed");
    }
    // Recover each face's coordinates directly from integer coordinates, using
    // a nonzero 2D minor; do not reuse producer barycentric propagation.
    for (U f = 0; f < nface; f++) {
      U root = f, depth = 0;
      while (faces[root].parent != NIL) {
        root = faces[root].parent;
        check(++depth <= 16, "cyclic or excessively deep physical hierarchy");
      }
      auto rv = faces[root].v;
      auto a = points[rv[0]].x, b = points[rv[1]].x, c = points[rv[2]].x;
      int j = -1, k = -1;
      I det = 0;
      for (int u = 0; u < 3; u++)
        for (int v = u + 1; v < 3; v++) {
          I d = (b[u] - a[u]) * (c[v] - a[v]) - (b[v] - a[v]) * (c[u] - a[u]);
          if (d) {
            j = u;
            k = v;
            det = d;
            break;
          }
        }
      check(j >= 0, "degenerate face");
      I scale = I(1) << depth;
      check(det % scale == 0, "barycentric determinant not divisible");
      I denom = det / scale;
      array<array<U, 3>, 3> bar;
      for (int v = 0; v < 3; v++) {
        auto x = points[faces[f].v[v]].x;
        I nb = (x[j] - a[j]) * (c[k] - a[k]) - (x[k] - a[k]) * (c[j] - a[j]);
        I nc = (b[j] - a[j]) * (x[k] - a[k]) - (b[k] - a[k]) * (x[j] - a[j]);
        check(nb % denom == 0 && nc % denom == 0,
              "non-dyadic barycentric coordinate");
        I ib = nb / denom, ic = nc / denom, ia = scale - ib - ic;
        check(ia >= 0 && ib >= 0 && ic >= 0 && ia <= scale && ib <= scale &&
                  ic <= scale,
              "face outside root triangle");
        bar[v] = {U(ia), U(ib), U(ic)};
        for (int h = 0; h < 3; h++)
          check(scale * x[h] == ia * a[h] + ib * b[h] + ic * c[h],
                "face off affine root plane");
      }
      U rep = root;
      vector<array<int, 3>> group{{0, 1, 2}};
      auto s = seeds.find(root);
      if (s != seeds.end()) {
        rep = s->second.rep;
        auto saved = bar;
        for (int v = 0; v < 3; v++)
          for (int h = 0; h < 3; h++)
            bar[v][h] = saved[v][s->second.frame[h]];
        group = s->second.group;
      }
      array<U, 9> best;
      best.fill(NIL);
      for (auto g : group) {
        array<array<U, 3>, 3> bb;
        for (int v = 0; v < 3; v++)
          for (int h = 0; h < 3; h++)
            bb[v][h] = bar[v][g[h]];
        sort(bb.begin(), bb.end());
        array<U, 9> flat;
        for (int v = 0; v < 3; v++)
          for (int h = 0; h < 3; h++)
            flat[3 * v + h] = bb[v][h];
        best = min(best, flat);
      }
      check(cells[faces[f].cell].key == Key{rep, depth, best},
            "canonical barycentric key mismatch");
      if (faces[f].child[0] != NIL) {
        array<U, 4> ch;
        for (int v = 0; v < 4; v++)
          ch[v] = faces[faces[f].child[v]].cell;
        sort(ch.begin(), ch.end());
        check(ch == cells[faces[f].cell].child, "canonical children mismatch");
      }
      if (f % 500000 == 0)
        cerr << "geometric faces " << f << '\n';
    }
    // Every canonical subdivision must have a physical source; forbid invented
    // descendants used solely to change prolongation rows.
    vector<uint8_t> source(ncell);
    for (auto f : faces)
      if (f.child[0] != NIL)
        source[f.cell] = 1;
    for (U c = 0; c < ncell; c++)
      check((cells[c].child[0] != NIL) == bool(source[c]),
            "canonical subdivision lacks physical witness");
    ifstream rows(prefix + ".rows.bin", ios::binary);
    uint64_t nrows;
    read(rows, nrows);
    check(nrows == 4 * nleaf, "row count mismatch");
    vector<uint8_t> identity(ndof);
    vector<U> incidence(nface);
    size_t nnz = 0, maxsupport = 0;
    uint64_t assembly_pairs = 0;
    for (size_t i = 0; i < nleaf; i++) {
      Tet t;
      array<U, 4> fs;
      U root, depth;
      uint64_t bits;
      read(topo, t);
      read(topo, fs);
      read(topo, root);
      read(topo, depth);
      read(topo, bits);
      check(root == plan[i].first && depth == plan[i].second.size(),
            "leaf plan mismatch");
      uint64_t expectedbits = 0;
      Tet expected = roots[root];
      for (char digit : plan[i].second) {
        int q = digit - '0';
        check(q >= 0 && q < 8, "invalid path digit");
        expectedbits = (expectedbits << 3) | q;
        U a = expected[0], b = expected[1], c = expected[2], d = expected[3];
        U ab = midpoint(a, b), ac = midpoint(a, c), ad = midpoint(a, d),
          bc = midpoint(b, c), bd = midpoint(b, d), cd = midpoint(c, d);
        array<Tet, 8> child = {Tet{a, ab, ac, ad},  Tet{b, ab, bc, bd},
                               Tet{c, ac, bc, cd},  Tet{d, ad, bd, cd},
                               Tet{ab, ac, ad, cd}, Tet{ab, ac, bc, cd},
                               Tet{ab, ad, bd, cd}, Tet{ab, bc, bd, cd}};
        expected = child[q];
      }
      check(expected == t && expectedbits == bits,
            "leaf differs from independent red refinement");
      set<U> localsupport;
      for (int omit = 0; omit < 4; omit++) {
        Tri tri;
        int z = 0;
        for (int k = 0; k < 4; k++)
          if (k != omit)
            tri[z++] = t[k];
        check(fs[omit] == lookupface(tri), "wrong local face");
        incidence[fs[omit]]++;
        map<U, Weight> expected;
        set<U> active;
        expectedrow(faces[fs[omit]].cell, 0, expected, active);
        U n;
        read(rows, n);
        check(n == expected.size(), "row support mismatch");
        for (auto kv : expected) {
          U dof, e;
          uint64_t num;
          read(rows, dof);
          read(rows, num);
          read(rows, e);
          check(dof == kv.first && Weight{num, e} == kv.second,
                "row weight mismatch");
          localsupport.insert(dof);
        }
        if (n == 1 && expected.begin()->second == Weight{1, 0})
          identity[expected.begin()->first] = 1;
        nnz += n;
        maxsupport = max(maxsupport, size_t(n));
      }
      assembly_pairs += uint64_t(localsupport.size()) * localsupport.size();
      if (i % 200000 == 0)
        cerr << "leaf replay " << i << '\n';
    }
    check(topo.peek() == EOF && rows.peek() == EOF,
          "unexpected trailing binary data");
    check(
        all_of(identity.begin(), identity.end(), [](auto x) { return x == 1; }),
        "prolongation rank witness missing");
    for (U f = 0; f < nface; f++)
      check(incidence[f] <= 2,
            "physical face incident to more than two leaves");
    ofstream report(prefix + ".verification.json");
    report << "{\n  \"schema\": \"d67-independent-face-moment-replay/v1\",\n  "
              "\"leaves\": "
           << nleaf << ",\n  \"master_dofs\": " << ndof
           << ",\n  \"prolongation_nonzeros\": " << nnz
           << ",\n  \"maximum_row_support\": " << maxsupport
           << ",\n  \"naive_local_assembly_pairs\": " << assembly_pairs
           << ",\n  \"geometric_barycentric_replay\": true,\n  "
              "\"leaf_path_replay\": true,\n  \"exact_prolongation_replay\": "
              "true,\n  \"full_column_rank_verified\": true,\n  "
              "\"source_geometry_binding_requires_wrapper\": true,\n  "
              "\"matrix_positivity_verified\": false,\n  "
              "\"spectral_exclusion_certified\": false\n}\n";
    cerr << "independent exact replay PASS\n";
    return 0;
  } catch (const exception &e) {
    cerr << e.what() << '\n';
    return 1;
  }
}
