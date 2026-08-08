#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

static int rank4(uint16_t m) {
  std::array<unsigned,4> rows{};
  for (int r=0;r<4;++r) rows[r]=(m>>(4*r))&15U;
  int rank=0;
  for (int b=0;b<4;++b) {
    int p=-1;
    for(int r=rank;r<4;++r) if((rows[r]>>b)&1U){p=r;break;}
    if(p<0) continue;
    std::swap(rows[rank],rows[p]);
    for(int r=0;r<4;++r) if(r!=rank && ((rows[r]>>b)&1U)) rows[r]^=rows[rank];
    ++rank;
  }
  return rank;
}

static uint16_t coset_min(uint16_t v, const std::vector<uint16_t>& span) {
  uint16_t z=v;
  for(auto s:span) z=std::min<uint16_t>(z,v^s);
  return z;
}

int main(int argc,char**argv){
  const char* out=(argc>1?argv[1]:"spread4.txt");
  std::array<uint8_t,65536> inv{};
  for(unsigned m=0;m<65536;++m) inv[m]=(rank4(m)==4);
  const uint16_t I=0x8421; // rows 0001,0010,0100,1000
  std::ofstream f(out);
  uint64_t count=0;
  std::vector<uint16_t> span{0,I};
  std::vector<uint16_t> gens{I};
  auto dfs = [&](auto&& self, uint16_t last)->void {
    if(gens.size()==4){
      ++count;
      f<<gens[0]<<' '<<gens[1]<<' '<<gens[2]<<' '<<gens[3]<<'\n';
      return;
    }
    for(unsigned vv=last+1;vv<65536;++vv){
      uint16_t v=(uint16_t)vv;
      if(coset_min(v,span)!=v) continue;
      bool ok=true;
      for(auto s:span) if(!inv[v^s]){ok=false;break;}
      if(!ok) continue;
      std::vector<uint16_t> old=span;
      span.reserve(old.size()*2);
      for(auto s:old) span.push_back(s^v);
      gens.push_back(v);
      self(self,v);
      gens.pop_back(); span=std::move(old);
    }
  };
  dfs(dfs,0);
  std::cerr<<"count "<<count<<"\n";
}
