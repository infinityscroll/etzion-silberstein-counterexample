#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

using U=uint16_t;
static U mul(U a,U b){U z=0;for(int i=0;i<4;++i)for(int j=0;j<4;++j){int x=0;for(int k=0;k<4;++k)x^=((a>>(4*i+k))&1)&((b>>(4*k+j))&1);z|=x<<(4*i+j);}return z;}
static U tr(U a){U z=0;for(int i=0;i<4;++i)for(int j=0;j<4;++j)z|=((a>>(4*i+j))&1)<<(4*j+i);return z;}
static U inverse(U a){std::array<unsigned,4>r{};for(int i=0;i<4;++i)r[i]=((a>>(4*i))&15)|(1U<<(4+i));for(int b=0;b<4;++b){int p=-1;for(int i=b;i<4;++i)if((r[i]>>b)&1){p=i;break;}if(p<0)return 0;std::swap(r[b],r[p]);for(int i=0;i<4;++i)if(i!=b&&((r[i]>>b)&1))r[i]^=r[b];}U z=0;for(int i=0;i<4;++i)z|=((r[i]>>4)&15)<<(4*i);return z;}
static int rank4(U m){std::array<unsigned,4>r{};for(int i=0;i<4;++i)r[i]=(m>>(4*i))&15;int k=0;for(int b=0;b<4;++b){int p=-1;for(int i=k;i<4;++i)if((r[i]>>b)&1){p=i;break;}if(p<0)continue;std::swap(r[k],r[p]);for(int i=k+1;i<4;++i)if((r[i]>>b)&1)r[i]^=r[k];++k;}return k;}

template<size_t N> static std::array<U,N> rref(std::vector<U> vs){int row=0;for(int b=15;b>=0;--b){int p=-1;for(int i=row;i<(int)vs.size();++i)if((vs[i]>>b)&1){p=i;break;}if(p<0)continue;std::swap(vs[row],vs[p]);for(int i=0;i<(int)vs.size();++i)if(i!=row&&((vs[i]>>b)&1))vs[i]^=vs[row];++row;}std::array<U,N>o{};int n=0;for(U v:vs)if(v)o[n++]=v;if(n!=(int)N){std::cerr<<"dimension "<<n<<" expected "<<N<<"\n";std::exit(2);}return o;}
template<size_t N> static std::vector<U> span(const std::array<U,N>&b){std::vector<U>w;w.reserve(1U<<N);for(unsigned x=0;x<(1U<<N);++x){U z=0;for(size_t i=0;i<N;++i)if((x>>i)&1)z^=b[i];w.push_back(z);}return w;}

struct DSU{std::vector<uint32_t>p,sz;explicit DSU(size_t n):p(n),sz(n,1){std::iota(p.begin(),p.end(),0);}uint32_t get(uint32_t x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}void add(uint32_t a,uint32_t b){a=get(a);b=get(b);if(a==b)return;if(sz[a]<sz[b])std::swap(a,b);p[b]=a;sz[a]+=sz[b];}};
static uint32_t span_key(unsigned x,unsigned y){std::array<unsigned,3>z{x,y,x^y};std::sort(z.begin(),z.end());return(z[0]<<12)|z[1];}

int main(int argc,char**argv){
  if(argc<3){std::cerr<<"usage: class(2|3) output.jsonl\n";return 1;}int cls=std::stoi(argv[1]);const char*out=argv[2];
  std::array<U,4>S = cls==2 ? std::array<U,4>{0x2ba7,0x8421,0xe395,0x16fb}:std::array<U,4>{0x263e,0x53b7,0x8421,0x3842};S=rref<4>(std::vector<U>(S.begin(),S.end()));auto sw=span(S);std::array<uint8_t,65536>ins{};for(U x:sw){ins[x]=1;if(x&&rank4(x)!=4){std::cerr<<"spread rank fail\n";return 3;}}
  std::vector<U>dw;for(unsigned x=0;x<65536;++x){bool ok=true;for(U s:S)if(std::popcount((unsigned)((U)x&s))&1){ok=false;break;}if(ok)dw.push_back((U)x);}auto D=rref<12>(dw);dw=span(D);std::array<int16_t,65536>msg;msg.fill(-1);for(unsigned x=0;x<4096;++x)msg[dw[x]]=x;
  std::array<unsigned,5>rd{};for(U x:dw)++rd[rank4(x)];std::cerr<<"class "<<cls<<" spread_basis";for(U x:S)std::cerr<<' '<<x;std::cerr<<" dual_basis";for(U x:D)std::cerr<<' '<<x;std::cerr<<" ranks";for(int i=0;i<5;++i)std::cerr<<' '<<rd[i];std::cerr<<"\n";if(rd[0]!=1||rd[1]!=0)return 4;

  std::vector<U>gl;std::array<U,65536>iv{};for(unsigned x=0;x<65536;++x)if(U y=inverse((U)x)){gl.push_back((U)x);iv[x]=y;}std::set<std::array<uint16_t,4096>> maps;
  for(U Q:gl){U qi=iv[Q];for(U N:sw)if(N){U P=mul(N,qi);bool ok=true;for(U s:sw)if(!ins[mul(mul(P,s),Q)]){ok=false;break;}if(!ok)continue;U pit=tr(iv[P]),qit=tr(qi);std::array<uint16_t,4096>T{};for(unsigned x=0;x<4096;++x){U y=mul(mul(pit,dw[x]),qit);if(msg[y]<0){std::cerr<<"dual aut fail\n";return 5;}T[x]=msg[y];}
      std::array<uint16_t,4096>invmsg{};for(unsigned x=0;x<4096;++x)invmsg[T[x]]=x;std::array<uint16_t,4096>fm{};for(unsigned l=0;l<4096;++l){unsigned lp=0;for(int j=0;j<12;++j)if(std::popcount(l&invmsg[1U<<j])&1)lp|=1U<<j;fm[l]=lp;}maps.insert(fm);
  }}
  std::cerr<<"automorphisms "<<maps.size()<<"\n";

  std::vector<uint16_t>aa,bb;std::vector<int32_t>ix(1U<<24,-1);aa.reserve(2794155);bb.reserve(2794155);for(unsigned a=1;a<4096;++a)for(unsigned b=a+1;b<4096;++b)if((a^b)>b){ix[(a<<12)|b]=aa.size();aa.push_back(a);bb.push_back(b);}DSU dsu(aa.size());
  for(const auto&fm:maps)for(uint32_t i=0;i<aa.size();++i){int32_t j=ix[span_key(fm[aa[i]],fm[bb[i]])];if(j<0)return 6;dsu.add(i,j);}
  std::vector<unsigned>r2;for(unsigned x=1;x<4096;++x)if(rank4(dw[x])==2)r2.push_back(x);std::vector<std::tuple<int,uint32_t,unsigned,unsigned,unsigned>>reps;uint64_t mass=0;for(uint32_t i=0;i<aa.size();++i)if(dsu.get(i)==i){int cnt=0;for(auto x:r2)if(!(std::popcount(x&aa[i])&1)&&!(std::popcount(x&bb[i])&1))++cnt;reps.emplace_back(cnt,dsu.sz[i],aa[i],bb[i],aa[i]^bb[i]);mass+=dsu.sz[i];}std::sort(reps.begin(),reps.end());std::ofstream f(out);for(auto[cnt,os,a,b,c]:reps)f<<"{\"rank2_count\":"<<cnt<<",\"orbit_size\":"<<os<<",\"span\":["<<a<<","<<b<<","<<c<<"]}\n";std::cerr<<"orbits "<<reps.size()<<" mass "<<mass<<"\n";
  int last=-1,n=0;uint64_t bm=0;for(auto[cnt,os,a,b,c]:reps){if(cnt!=last){if(last>=0)std::cerr<<"count "<<last<<" orbits "<<n<<" mass "<<bm<<"\n";last=cnt;n=0;bm=0;}++n;bm+=os;}if(last>=0)std::cerr<<"count "<<last<<" orbits "<<n<<" mass "<<bm<<"\n";
}
