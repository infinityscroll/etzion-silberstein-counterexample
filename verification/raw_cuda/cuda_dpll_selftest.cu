#include <cuda_runtime.h>

#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

static constexpr int MAX_CLAUSES=128;

__device__ __forceinline__ void reduce_form(uint64_t&m,int&r,const uint64_t*eq,uint64_t rb){while(m){int p=63-__clzll(m);if(!eq[p])break;m^=eq[p];r^=(rb>>p)&1ULL;}}
__device__ __forceinline__ bool add_eq(uint64_t m,int r,uint64_t*eq,uint64_t&rb,uint8_t*trail,int&nt){reduce_form(m,r,eq,rb);if(!m)return r==0;int p=63-__clzll(m);if(nt>=48)return false;eq[p]=m;if(r)rb|=1ULL<<p;else rb&=~(1ULL<<p);trail[nt++]=p;return true;}
__device__ __forceinline__ void restore(int mark,uint64_t*eq,uint64_t&rb,uint8_t*trail,int&nt){while(nt>mark){int p=trail[--nt];eq[p]=0;rb&=~(1ULL<<p);}}

__device__ uint8_t solve(const uint64_t*ca,const uint64_t*cb,int nc,uint64_t&nodes){
  uint64_t eq[48]={},rb=0,fa[49],fb[49];uint8_t trail[48],fva[49],fvb[49],fm[49],fs[49];int nt=0,depth=0;nodes=0;
  while(true){++nodes;bool conflict=false,have=false;uint64_t xa=0,xb=0;int va=0,vb=0,bn=99,bw=999;
    for(int i=0;i<nc;++i){uint64_t a=ca[i],b=cb[i];int ar=0,br=0;reduce_form(a,ar,eq,rb);reduce_form(b,br,eq,rb);if((!a&&ar)||(!b&&br))continue;if(!a&&!ar&&!b&&!br){conflict=true;break;}int n=(a!=0)+(b!=0),w=__popcll(a)+__popcll(b);if(!have||n<bn||(n==bn&&w<bw)){have=true;bn=n;bw=w;xa=a;xb=b;va=ar;vb=br;}}
    if(!conflict&&!have)return 1;
    if(!conflict){if(depth>=49)return 3;fa[depth]=xa;fb[depth]=xb;fva[depth]=va;fvb[depth]=vb;fm[depth]=nt;if(!xa){fs[depth++]=2;if(add_eq(xb,1^vb,eq,rb,trail,nt))continue;}else if(!xb){fs[depth++]=2;if(add_eq(xa,1^va,eq,rb,trail,nt))continue;}else{fs[depth++]=1;if(add_eq(xa,1^va,eq,rb,trail,nt))continue;}}
    bool down=false;while(depth>0){int f=depth-1;restore(fm[f],eq,rb,trail,nt);if(fs[f]==1){fs[f]=2;bool ok=add_eq(fa[f],fva[f],eq,rb,trail,nt);if(ok)ok=add_eq(fb[f],1^fvb[f],eq,rb,trail,nt);if(ok){down=true;break;}restore(fm[f],eq,rb,trail,nt);}--depth;}if(down)continue;return 0;
  }
}

__global__ void test_kernel(const uint64_t*ca,const uint64_t*cb,const uint16_t*nc,int tests,uint8_t*status,uint64_t*nodes){int t=blockIdx.x*blockDim.x+threadIdx.x;if(t>=tests)return;status[t]=solve(ca+t*MAX_CLAUSES,cb+t*MAX_CLAUSES,nc[t],nodes[t]);}

static bool truth(const uint64_t*ca,const uint64_t*cb,int nc,int variables){for(uint64_t x=0;x<(1ULL<<variables);++x){bool ok=true;for(int i=0;i<nc;++i)if(!(std::popcount(x&ca[i])&1)&&!(std::popcount(x&cb[i])&1)){ok=false;break;}if(ok)return true;}return false;}
static void check(cudaError_t e,const char*w){if(e!=cudaSuccess){std::cerr<<w<<": "<<cudaGetErrorString(e)<<"\n";std::exit(2);}}

int main(int argc,char**argv){std::string out=argc>1?argv[1]:"cuda_dpll_selftest.json";constexpr int tests=256;constexpr uint64_t seed=8808001;std::mt19937_64 rng(seed);std::vector<uint64_t>ca(tests*MAX_CLAUSES),cb(tests*MAX_CLAUSES);std::vector<uint16_t>nc(tests);std::vector<uint8_t>expected(tests);
  for(int t=0;t<tests;++t){int n=4+rng()%11;int m=(t%2==0)?std::max(1,n/2):std::min(MAX_CLAUSES,8*n);nc[t]=m;uint64_t lim=(1ULL<<n)-1;for(int i=0;i<m;++i){uint64_t a=rng()&lim,b=rng()&lim;if(!a)a=1;if(!b)b=1;ca[t*MAX_CLAUSES+i]=a;cb[t*MAX_CLAUSES+i]=b;}expected[t]=truth(ca.data()+t*MAX_CLAUSES,cb.data()+t*MAX_CLAUSES,m,n);}
  uint64_t *da,*db,*dnodes;uint16_t*dnc;uint8_t*ds;check(cudaMalloc(&da,ca.size()*8),"malloc a");check(cudaMalloc(&db,cb.size()*8),"malloc b");check(cudaMalloc(&dnc,nc.size()*2),"malloc nc");check(cudaMalloc(&ds,tests),"malloc status");check(cudaMalloc(&dnodes,tests*8),"malloc nodes");check(cudaMemcpy(da,ca.data(),ca.size()*8,cudaMemcpyHostToDevice),"copy a");check(cudaMemcpy(db,cb.data(),cb.size()*8,cudaMemcpyHostToDevice),"copy b");check(cudaMemcpy(dnc,nc.data(),nc.size()*2,cudaMemcpyHostToDevice),"copy nc");test_kernel<<<(tests+127)/128,128>>>(da,db,dnc,tests,ds,dnodes);check(cudaDeviceSynchronize(),"sync");std::vector<uint8_t>got(tests);std::vector<uint64_t>nodes(tests);check(cudaMemcpy(got.data(),ds,tests,cudaMemcpyDeviceToHost),"get status");check(cudaMemcpy(nodes.data(),dnodes,tests*8,cudaMemcpyDeviceToHost),"get nodes");int mismatch=0,sat=0,unsat=0,errors=0;uint64_t total_nodes=0;for(int i=0;i<tests;++i){if(got[i]>1)++errors;else if(got[i])++sat;else ++unsat;if(got[i]!=(uint8_t)expected[i])++mismatch;total_nodes+=nodes[i];}std::ofstream f(out);f<<"{\n  \"seed\": "<<seed<<",\n  \"tests\": "<<tests<<",\n  \"variable_range\": [4, 14],\n  \"SAT\": "<<sat<<",\n  \"UNSAT\": "<<unsat<<",\n  \"errors\": "<<errors<<",\n  \"truth_table_mismatches\": "<<mismatch<<",\n  \"total_gpu_nodes\": "<<total_nodes<<"\n}\n";std::cout<<"tests "<<tests<<" SAT "<<sat<<" UNSAT "<<unsat<<" errors "<<errors<<" mismatches "<<mismatch<<" nodes "<<total_nodes<<"\n";return mismatch||errors?1:0;}
