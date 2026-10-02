// Exact adaptive CR face hierarchy. This is a producer, not a spectral verifier.
// Build: g++ -std=c++17 -O2 moment_map.cpp -o /tmp/d67-moment-map
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using I=__int128_t; using U=uint32_t;
using Tri=array<U,3>; using Tet=array<U,4>;
const U NIL=UINT32_MAX;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
I parse(const string&s){I x=0;for(char c:s)if(c!='-')x=10*x+c-'0';return s[0]=='-'?-x:x;}
struct Node {int64_t a,b;I r;bool operator==(const Node&v)const{return a==v.a&&b==v.b&&r==v.r;}};
uint64_t mix(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
struct NH{size_t operator()(const Node&n)const{return mix(n.a)^mix(n.b+7)^mix(uint64_t(n.r))^mix(uint64_t(n.r>>64)+13);}};
struct TH{size_t operator()(const Tri&t)const{return mix(t[0])^mix(uint64_t(t[1])+10000000000ULL)^mix(uint64_t(t[2])+20000000000ULL);}};
Tri sorttri(Tri t){sort(t.begin(),t.end());return t;}
struct Face {Tri v;U parent=NIL;array<U,4>child={NIL,NIL,NIL,NIL};U cell=NIL;};
struct Seed {U rep;array<int,3>frame;vector<array<int,3>>group;};
// Canonical coordinates are sorted barycentric rows on a representative root.
struct Key {U root;U depth;array<U,9>b;bool operator==(const Key&k)const{return root==k.root&&depth==k.depth&&b==k.b;}};
struct KH {size_t operator()(const Key&k)const{uint64_t h=mix(k.root)^mix(k.depth+991);for(auto v:k.b)h=mix(h^v);return h;}};
struct Cell {Key key;array<U,4>child={NIL,NIL,NIL,NIL};U dof=NIL;};
vector<Node>nodes; vector<array<U,2>>nodeparents; unordered_map<Node,U,NH>nodeindex;
vector<Face>faces;unordered_map<Tri,U,TH>faceindex;
vector<Cell>cells;unordered_map<Key,U,KH>cellindex;
unordered_map<U,Seed>seeds;
vector<Tet>roots,leaves;vector<array<U,4>>leaffaces;vector<U>leafroots;vector<uint64_t>leafpaths;vector<U>leafdepth;
vector<pair<U,string>>plan;size_t pos=0;
U face(Tri v){v=sorttri(v);auto p=faceindex.emplace(v,faces.size());if(p.second)faces.push_back({v});return p.first->second;}
array<U,4> tetfaces(const Tet&t){array<U,4>f;for(int i=0;i<4;i++){Tri v;int k=0;for(int j=0;j<4;j++)if(i!=j)v[k++]=t[j];f[i]=face(v);}return f;}
U mid(U x,U y){if(x>y)swap(x,y);auto a=nodes[x],b=nodes[y];I ar=I(a.a)+b.a,br=I(a.b)+b.b,rr=a.r+b.r;need(ar%2==0&&br%2==0&&rr%2==0,"midpoint outside grid");Node n{int64_t(ar/2),int64_t(br/2),rr/2};auto p=nodeindex.emplace(n,nodes.size());if(p.second){nodes.push_back(n);nodeparents.push_back({x,y});}return p.first->second;}
void splitface(U id){Tri v=faces[id].v;U x=v[0],y=v[1],z=v[2],xy=mid(x,y),xz=mid(x,z),yz=mid(y,z);array<Tri,4> ts={Tri{x,xy,xz},Tri{y,xy,yz},Tri{z,xz,yz},Tri{xy,xz,yz}};array<U,4>ch;for(int i=0;i<4;i++){ch[i]=face(ts[i]);need(faces[ch[i]].parent==NIL||faces[ch[i]].parent==id,"inconsistent face ancestry");faces[ch[i]].parent=id;}if(faces[id].child[0]!=NIL)need(faces[id].child==ch,"inconsistent face subdivision");faces[id].child=ch;}
void walk(U root,const Tet&t,const string&path,uint64_t bits){need(pos<plan.size()&&plan[pos].first==root,"incomplete root tree");if(plan[pos].second==path){leaves.push_back(t);leaffaces.push_back(tetfaces(t));leafroots.push_back(root);leafpaths.push_back(bits);leafdepth.push_back(path.size());pos++;return;}need(plan[pos].second.compare(0,path.size(),path)==0&&path.size()<16,"invalid prefix tree");auto fs=tetfaces(t);for(U f:fs)splitface(f);U a=t[0],b=t[1],c=t[2],d=t[3];U ab=mid(a,b),ac=mid(a,c),ad=mid(a,d),bc=mid(b,c),bd=mid(b,d),cd=mid(c,d);array<Tet,8> ch={Tet{a,ab,ac,ad},Tet{b,ab,bc,bd},Tet{c,ac,bc,cd},Tet{d,ad,bd,cd},Tet{ab,ac,ad,cd},Tet{ab,ac,bc,cd},Tet{ab,ad,bd,cd},Tet{ab,bc,bd,cd}};for(auto&v:ch)tetfaces(v);for(int k=0;k<8;k++)walk(root,ch[k],path+char('0'+k),(bits<<3)|k);}
// Recursion supplies barycentric coordinates; a separate verifier must recover
// them geometrically, so it does not rely on this producer's recurrence.
struct Frame {U root,depth;array<array<U,3>,3>b;};
vector<Frame>frames;
Frame getframe(U id){if(frames[id].root!=NIL)return frames[id];auto f=faces[id];Frame out;out.depth=0;out.root=id;out.b={array<U,3>{1,0,0},array<U,3>{0,1,0},array<U,3>{0,0,1}};if(f.parent!=NIL){auto p=getframe(f.parent);out.root=p.root;out.depth=p.depth+1;auto pv=faces[f.parent].v;for(int i=0;i<3;i++){U n=f.v[i];int v=-1;for(int j=0;j<3;j++)if(n==pv[j])v=j;if(v>=0){for(int j=0;j<3;j++)out.b[i][j]=2*p.b[v][j];}else{bool found=false;for(int j=0;j<3;j++)for(int k=j+1;k<3;k++)if(mid(pv[j],pv[k])==n){for(int l=0;l<3;l++)out.b[i][l]=p.b[j][l]+p.b[k][l];found=true;}need(found,"child vertex not a face midpoint");}}}frames[id]=out;return out;}
U canonical(U f){auto fr=getframe(f);Key k;k.root=fr.root;k.depth=fr.depth;array<int,3>identity={0,1,2};vector<array<int,3>>group{identity};auto b=fr.b;auto s=seeds.find(fr.root);if(s!=seeds.end()){k.root=s->second.rep;for(int i=0;i<3;i++)for(int j=0;j<3;j++)b[i][j]=fr.b[i][s->second.frame[j]];group=s->second.group;}array<U,9>best;best.fill(UINT32_MAX);for(auto g:group){array<array<U,3>,3>tmp;for(int i=0;i<3;i++)for(int j=0;j<3;j++)tmp[i][j]=b[i][g[j]];sort(tmp.begin(),tmp.end());array<U,9>flat;for(int i=0;i<3;i++)for(int j=0;j<3;j++)flat[3*i+j]=tmp[i][j];best=min(best,flat);}k.b=best;auto p=cellindex.emplace(k,cells.size());if(p.second)cells.push_back({k});return p.first->second;}
using Weight=pair<uint64_t,U>;
void addweight(Weight&to,Weight from){U exp=max(to.second,from.second);need(exp<=60,"weight denominator overflow");uint64_t a=to.first<<(exp-to.second),b=from.first<<(exp-from.second);need(UINT64_MAX-a>=b,"weight numerator overflow");to={a+b,exp};while(to.second&&!(to.first&1)){to.first>>=1;to.second--;}}
void expand(U c,U exp,map<U,Weight>&w){auto&v=cells[c];if(v.child[0]==NIL){addweight(w[v.dof],{1,exp});return;}for(U child:v.child)expand(child,exp+2,w);}
template<class T>void put(ofstream&f,const T&x){f.write(reinterpret_cast<const char*>(&x),sizeof(x));}
int main(int argc,char**argv){try{need(argc==3,"usage: moment_map INPUT OUTPUT_PREFIX");ifstream in(argv[1]);need(bool(in),"input open failed");size_t nn,nr,ns,nl;string zs,rs;in>>nn>>nr>>ns>>nl>>zs>>rs;nodes.reserve(nl);nodeparents.reserve(nl);nodeindex.reserve(nl);faces.reserve(4*nl);faceindex.reserve(4*nl);for(size_t i=0;i<nn;i++){int64_t a,b;string r;in>>a>>b>>r;Node n{a,b,parse(r)};need(nodeindex.emplace(n,nodes.size()).second,"duplicate initial node");nodes.push_back(n);nodeparents.push_back({NIL,NIL});}roots.resize(nr);for(auto&t:roots){for(auto&n:t)in>>n;tetfaces(t);}for(size_t i=0;i<ns;i++){Tri src,rep;for(auto&n:src)in>>n;for(auto&n:rep)in>>n;Seed s;s.rep=face(rep);for(auto&n:s.frame)in>>n;size_t ng;in>>ng;s.group.resize(ng);for(auto&g:s.group)for(auto&n:g)in>>n;need(seeds.emplace(face(src),s).second,"duplicate root seed");}plan.resize(nl);for(auto&p:plan){in>>p.first>>p.second;if(p.second=="*")p.second="";}need(bool(in),"truncated input");for(U r=0;r<nr;r++){walk(r,roots[r],"",0);if(r%500==0)cerr<<"roots "<<r<<" faces "<<faces.size()<<'\n';}need(pos==plan.size(),"unused leaves");frames.resize(faces.size());for(auto&f:frames)f.root=NIL;cells.reserve(faces.size());cellindex.reserve(faces.size());for(U f=0;f<faces.size();f++)faces[f].cell=canonical(f);need(frames.size()==faces.size(),"canonicalization created extra face");for(U f=0;f<faces.size();f++)if(faces[f].child[0]!=NIL){array<U,4>ch;for(int j=0;j<4;j++)ch[j]=faces[faces[f].child[j]].cell;sort(ch.begin(),ch.end());auto&c=cells[faces[f].cell];if(c.child[0]!=NIL)need(c.child==ch,"paired child partitions differ");c.child=ch;}
U ndof=0;for(auto&c:cells)if(c.child[0]==NIL)c.dof=ndof++;vector<uint8_t>witness(ndof);vector<U>incidence(ndof);size_t nz=0,maxrow=0;string prefix=argv[2];ofstream rows(prefix+".rows.bin",ios::binary);uint64_t nrows=4*leaves.size();put(rows,nrows);for(auto fs:leaffaces)for(U f:fs){map<U,Weight>w;expand(faces[f].cell,0,w);Weight sum{0,0};for(auto kv:w)addweight(sum,kv.second);need(sum.first==1&&sum.second==0,"row weights fail to sum to one");if(w.size()==1&&w.begin()->second==Weight{1,0})witness[w.begin()->first]=1;for(auto kv:w)incidence[kv.first]++;U sz=w.size();put(rows,sz);for(auto kv:w){put(rows,kv.first);put(rows,kv.second.first);put(rows,kv.second.second);}nz+=w.size();maxrow=max(maxrow,w.size());}need(all_of(witness.begin(),witness.end(),[](auto x){return x==1;}),"master lacks identity-row rank witness");rows.close();
ofstream topo(prefix+".topology.bin",ios::binary);uint64_t magic=0x4436374d4f4d3031ULL;put(topo,magic);for(uint64_t n:{nodes.size(),faces.size(),cells.size(),leaves.size(),size_t(ndof)})put(topo,n);for(size_t i=0;i<nodes.size();i++){put(topo,nodes[i].a);put(topo,nodes[i].b);uint64_t lo=uint64_t(nodes[i].r),hi=uint64_t(nodes[i].r>>64);put(topo,lo);put(topo,hi);put(topo,nodeparents[i]);}for(auto f:faces){put(topo,f.v);put(topo,f.parent);put(topo,f.child);put(topo,f.cell);}for(auto c:cells){put(topo,c.key.root);put(topo,c.key.depth);put(topo,c.key.b);put(topo,c.child);put(topo,c.dof);}for(size_t i=0;i<leaves.size();i++){put(topo,leaves[i]);put(topo,leaffaces[i]);put(topo,leafroots[i]);put(topo,leafdepth[i]);put(topo,leafpaths[i]);}topo.close();ofstream summary(prefix+".json");summary<<"{\n  \"schema\": \"d67-face-moment-producer/v1\",\n  \"nodes\": "<<nodes.size()<<",\n  \"physical_faces\": "<<faces.size()<<",\n  \"canonical_cells\": "<<cells.size()<<",\n  \"leaves\": "<<leaves.size()<<",\n  \"master_dofs\": "<<ndof<<",\n  \"prolongation_nonzeros\": "<<nz<<",\n  \"maximum_row_support\": "<<maxrow<<",\n  \"all_rows_sum_to_one\": true,\n  \"all_master_identity_witnesses\": true,\n  \"independently_verified\": false,\n  \"spectral_exclusion_certified\": false\n}\n";cerr<<"complete: "<<ndof<<" dofs, max row "<<maxrow<<'\n';return 0;}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
