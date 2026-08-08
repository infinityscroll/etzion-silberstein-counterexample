#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <vector>

static unsigned gf_mul(unsigned a,unsigned b){unsigned r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&16)a=(a&15)^3;}return r;}
static unsigned matrix_of(unsigned x){unsigned a=x&15,b=(x>>4)&15,c=(x>>8)&15,m=0;for(int j=0;j<4;++j){unsigned z=1U<<j,z2=gf_mul(z,z),z4=gf_mul(z2,z2);m|=(gf_mul(a,z)^gf_mul(b,z2)^gf_mul(c,z4))<<(4*j);}return m;}
static std::array<unsigned,4> rows_of(unsigned m){std::array<unsigned,4> r{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)r[i]|=((m>>(4*j+i))&1U)<<j;return r;}
static int rank4(unsigned m){auto r=rows_of(m);int k=0;for(int b=0;b<4;++b){int p=-1;for(int i=k;i<4;++i)if((r[i]>>b)&1U){p=i;break;}if(p<0)continue;std::swap(r[k],r[p]);for(int i=k+1;i<4;++i)if((r[i]>>b)&1U)r[i]^=r[k];++k;}return k;}

static std::vector<unsigned> vector_basis(const std::vector<unsigned>& vs){
  std::array<unsigned,12> p{};
  for(unsigned v:vs){while(v){int b=31-std::countl_zero(v);if(p[b])v^=p[b];else{p[b]=v;break;}}}
  std::vector<unsigned> out;for(int b=11;b>=0;--b)if(p[b])out.push_back(p[b]);return out;
}

struct EquationState{std::array<uint64_t,40> m{};std::array<uint8_t,40> r{};};
static std::pair<uint64_t,int> reduce(uint64_t m,int r,const EquationState& e){while(m){int p=63-std::countl_zero(m);if(!e.m[p])break;m^=e.m[p];r^=e.r[p];}return {m,r};}
static bool add(EquationState& e,uint64_t m,int r){auto z=reduce(m,r,e);m=z.first;r=z.second;if(!m)return r==0;int p=63-std::countl_zero(m);e.m[p]=m;e.r[p]=(uint8_t)r;return true;}

struct Solver{
  std::vector<std::pair<uint64_t,uint64_t>> cs;uint64_t nodes=0,conflicts=0,answer=0;
  bool dfs(const EquationState& e){
    ++nodes;bool have=false;uint64_t ca=0,cb=0;int va=0,vb=0,bestn=99,bestw=999;
    for(auto [a,b]:cs){auto ra=reduce(a,0,e),rb=reduce(b,0,e);
      if((!ra.first&&ra.second)||(!rb.first&&rb.second))continue;
      if(!ra.first&&!ra.second&&!rb.first&&!rb.second){++conflicts;return false;}
      int n=(ra.first!=0)+(rb.first!=0),w=std::popcount(ra.first)+std::popcount(rb.first);
      if(!have||std::pair(n,w)<std::pair(bestn,bestw)){have=true;bestn=n;bestw=w;ca=ra.first;va=ra.second;cb=rb.first;vb=rb.second;}
    }
    if(!have){uint64_t x=0;for(int p=0;p<40;++p)if(e.m[p]){uint64_t low=e.m[p]&((1ULL<<p)-1);int v=e.r[p]^(std::popcount(x&low)&1);if(v)x|=1ULL<<p;}answer=x;return true;}
    if(!ca){EquationState q=e;if(!add(q,cb,1^vb))return false;return dfs(q);}
    if(!cb){EquationState q=e;if(!add(q,ca,1^va))return false;return dfs(q);}
    EquationState l=e;if(add(l,ca,1^va)&&dfs(l))return true;
    EquationState r=e;if(add(r,ca,va)&&add(r,cb,1^vb)&&dfs(r))return true;
    return false;
  }
};

static std::array<unsigned,2> null_basis(unsigned matrix){
  auto rows=rows_of(matrix);std::vector<unsigned> z;
  for(unsigned v=1;v<16;++v){bool ok=true;for(auto r:rows)if(std::popcount(r&v)&1){ok=false;break;}if(ok)z.push_back(v);}
  auto b=vector_basis(z);if(b.size()!=2){std::cerr<<"bad null\n";std::exit(3);}return {b[0],b[1]};
}

struct Instance{std::vector<unsigned> basis;std::vector<std::pair<uint64_t,uint64_t>> cs;};
static Instance make_instance(unsigned l1,unsigned l2,const std::array<unsigned,4096>& mats,const std::array<uint8_t,4096>& ranks){
  std::vector<unsigned> words;words.reserve(1024);for(unsigned x=0;x<4096;++x)if(!(std::popcount(x&l1)&1)&&!(std::popcount(x&l2)&1))words.push_back(x);
  auto basis=vector_basis(words);if(basis.size()!=10){std::cerr<<"bad kernel\n";std::exit(4);}std::array<unsigned,4096> coord{};
  for(unsigned u=0;u<1024;++u){unsigned x=0;for(int i=0;i<10;++i)if((u>>i)&1)x^=basis[i];coord[x]=u;}
  std::set<std::pair<uint64_t,uint64_t>> uniq;
  for(unsigned x:words)if(x&&ranks[x]==2){auto ns=null_basis(mats[x]);uint64_t a=0,b=0;for(int j=0;j<4;++j){if((ns[0]>>j)&1)a|=(uint64_t)coord[x]<<(10*j);if((ns[1]>>j)&1)b|=(uint64_t)coord[x]<<(10*j);}if(a>b)std::swap(a,b);uniq.emplace(a,b);}
  return {basis,{uniq.begin(),uniq.end()}};
}

int main(int argc,char**argv){
  if(argc<6){std::cerr<<"usage: reps.jsonl lane lanes max_count out.jsonl [hit.json]\n";return 1;}
  std::string repfile=argv[1],outfile=argv[5],hitfile=argc>6?argv[6]:"orbit_hit.json";int lane=std::stoi(argv[2]),lanes=std::stoi(argv[3]),maxcnt=std::stoi(argv[4]);
  std::array<unsigned,4096> mats{};std::array<uint8_t,4096> ranks{};for(unsigned x=0;x<4096;++x){mats[x]=matrix_of(x);ranks[x]=rank4(mats[x]);}
  std::ifstream f(repfile);std::ofstream out(outfile);std::string s;std::regex re("rank2_count\\\":([0-9]+).*orbit_size\\\":([0-9]+).*span\\\":\\[([0-9]+),([0-9]+),([0-9]+)\\]");std::smatch m;unsigned idx=0,done=0;
  while(std::getline(f,s))if(std::regex_search(s,m,re)){int cnt=std::stoi(m[1]);unsigned os=std::stoul(m[2]),l1=std::stoul(m[3]),l2=std::stoul(m[4]),l3=std::stoul(m[5]);if(cnt>maxcnt)break;unsigned cur=idx++;if(cur%lanes!=(unsigned)lane)continue;
    auto t0=std::chrono::steady_clock::now();auto in=make_instance(l1,l2,mats,ranks);Solver sol;sol.cs=std::move(in.cs);EquationState e;bool sat=sol.dfs(e);double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    out<<"{\"rep_index\":"<<cur<<",\"rank2_count\":"<<cnt<<",\"orbit_size\":"<<os<<",\"span\":["<<l1<<","<<l2<<","<<l3<<"],\"constraints\":"<<sol.cs.size()<<",\"status\":\""<<(sat?"SAT":"UNSAT")<<"\",\"nodes\":"<<sol.nodes<<",\"conflicts\":"<<sol.conflicts<<",\"elapsed_s\":"<<sec;
    if(sat)out<<",\"lift_assignment\":"<<sol.answer;
    out<<"}\n";out.flush();++done;
    if(sat){std::ofstream h(hitfile);h<<"{\"span\":["<<l1<<","<<l2<<","<<l3<<"],\"lift_assignment\":"<<sol.answer<<"}\n";std::cerr<<"SAT rep "<<cur<<"\n";return 0;}
  }
  std::cerr<<"done "<<done<<"\n";
}
