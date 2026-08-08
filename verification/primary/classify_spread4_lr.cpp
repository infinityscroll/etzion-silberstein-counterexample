#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>

using U=uint16_t;
static U mul(U a,U b){U z=0;for(int i=0;i<4;++i)for(int j=0;j<4;++j){int x=0;for(int k=0;k<4;++k)x^=((a>>(4*i+k))&1)&((b>>(4*k+j))&1);z|=x<<(4*i+j);}return z;}
static int rank4(U m){std::array<unsigned,4>r{};for(int i=0;i<4;++i)r[i]=(m>>(4*i))&15;int k=0;for(int b=0;b<4;++b){int p=-1;for(int i=k;i<4;++i)if((r[i]>>b)&1){p=i;break;}if(p<0)continue;std::swap(r[k],r[p]);for(int i=0;i<4;++i)if(i!=k&&((r[i]>>b)&1))r[i]^=r[k];++k;}return k;}
static U inverse(U a){std::array<unsigned,4> r{};for(int i=0;i<4;++i)r[i]=((a>>(4*i))&15)|(1U<<(4+i));for(int b=0;b<4;++b){int p=-1;for(int i=b;i<4;++i)if((r[i]>>b)&1){p=i;break;}if(p<0)return 0;std::swap(r[b],r[p]);for(int i=0;i<4;++i)if(i!=b&&((r[i]>>b)&1))r[i]^=r[b];}U z=0;for(int i=0;i<4;++i)z|=((r[i]>>4)&15)<<(4*i);return z;}
static std::array<U,65536> make_inv(){std::array<U,65536> tab{};for(unsigned a=0;a<65536;++a)tab[a]=inverse((U)a);return tab;}
static U trans(U a){U z=0;for(int i=0;i<4;++i)for(int j=0;j<4;++j)z|=((a>>(4*i+j))&1)<<(4*j+i);return z;}

static std::array<U,4> basis(std::vector<U> vs){
  int row=0;for(int b=15;b>=0;--b){int p=-1;for(int i=row;i<(int)vs.size();++i)if((vs[i]>>b)&1){p=i;break;}if(p<0)continue;std::swap(vs[row],vs[p]);for(int i=0;i<(int)vs.size();++i)if(i!=row&&((vs[i]>>b)&1))vs[i]^=vs[row];++row;}
  std::array<U,4> out{};int n=0;for(U v:vs)if(v)out[n++]=v;if(n!=4){std::cerr<<"dim "<<n<<"\n";std::exit(3);}return out;
}
static uint64_t key(const std::array<U,4>&b){uint64_t k=0;for(int i=0;i<4;++i)k|=(uint64_t)b[i]<<(16*i);return k;}
static std::vector<U> words(const std::array<U,4>&b){std::vector<U>w;for(int x=0;x<16;++x){U z=0;for(int i=0;i<4;++i)if((x>>i)&1)z^=b[i];w.push_back(z);}return w;}
struct DSU{std::vector<int>p,s;DSU(int n):p(n),s(n,1){std::iota(p.begin(),p.end(),0);}int get(int x){return p[x]==x?x:p[x]=get(p[x]);}void add(int a,int b){a=get(a);b=get(b);if(a==b)return;if(s[a]<s[b])std::swap(a,b);p[b]=a;s[a]+=s[b];}};

int main(int argc,char**argv){const char* fn=argc>1?argv[1]:"spread4.txt";auto invtab=make_inv();std::ifstream f(fn);std::vector<std::array<U,4>> spaces;unsigned a,b,c,d;while(f>>a>>b>>c>>d)spaces.push_back(basis({(U)a,(U)b,(U)c,(U)d}));std::cerr<<"spaces "<<spaces.size()<<"\n";
  std::unordered_map<uint64_t,int> ix;for(int i=0;i<(int)spaces.size();++i)ix[key(spaces[i])]=i;if(ix.size()!=spaces.size()){std::cerr<<"duplicates\n";return 4;}DSU ds(spaces.size());
  std::vector<U> qs;for(int i=0;i<4;++i)for(int j=0;j<4;++j)if(i!=j)qs.push_back((U)(0x8421^(1U<<(4*i+j)))); // transvections
  for(int i=0;i<(int)spaces.size();++i){auto ws=words(spaces[i]);
    for(U m:ws)if(m){U mi=invtab[m];std::vector<U> z;for(U x:ws)z.push_back(mul(mi,x));auto it=ix.find(key(basis(z)));if(it==ix.end()){std::cerr<<"missing left\n";return 5;}ds.add(i,it->second);}
    for(U q:qs){U qi=invtab[q];std::vector<U>z;for(U x:ws)z.push_back(mul(mul(qi,x),q));auto it=ix.find(key(basis(z)));if(it==ix.end()){std::cerr<<"missing conj\n";return 6;}ds.add(i,it->second);}
  }
  std::vector<int> reps;for(int i=0;i<(int)spaces.size();++i)if(ds.get(i)==i)reps.push_back(i);std::cout<<"lr_orbits "<<reps.size()<<"\n";for(int r:reps){std::cout<<"rep "<<r<<" normalized_spaces "<<ds.s[r]<<" basis";for(U x:spaces[r])std::cout<<' '<<x;std::cout<<"\n";}
  std::cout<<"transpose_map\n";for(int r:reps){std::vector<U>z;for(U x:words(spaces[r]))z.push_back(trans(x));int j=ix.at(key(basis(z)));std::cout<<r<<" -> "<<ds.get(j)<<"\n";}
  std::cout<<"named_representatives\n";
  const std::array<std::pair<const char*,std::array<U,4>>,3> named{{
    {"field_dual_of_standard_D",{33439,16862,10743,7658}},
    {"published_class2",{0x2ba7,0x8421,0xe395,0x16fb}},
    {"published_class3",{0x263e,0x53b7,0x8421,0x3842}}
  }};
  for(const auto& [name,g]:named){auto bg=basis(std::vector<U>(g.begin(),g.end()));auto wg=words(bg);U m=*std::find_if(wg.begin(),wg.end(),[](U x){return x!=0;});std::vector<U>ng;for(U x:wg)ng.push_back(mul(invtab[m],x));int j=ix.at(key(basis(ng)));std::cout<<name<<" normalized_index "<<j<<" root "<<ds.get(j)<<" root_mass "<<ds.s[ds.get(j)]<<"\n";}
}
