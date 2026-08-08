#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Pair { uint16_t a, b; };
struct RankTwo { uint16_t message; uint8_t n1, n2; };

__device__ __forceinline__ void reduce_form(uint64_t &mask, int &rhs,
                                             const uint64_t *eq,
                                             uint64_t rhs_bits) {
  while (mask) {
    int pivot = 63 - __clzll(mask);
    uint64_t row = eq[pivot];
    if (!row) break;
    mask ^= row;
    rhs ^= int((rhs_bits >> pivot) & 1ULL);
  }
}

__device__ __forceinline__ bool add_equation(uint64_t mask, int rhs,
                                              uint64_t *eq,
                                              uint64_t &rhs_bits,
                                              uint8_t *trail,
                                              int &trail_size) {
  reduce_form(mask, rhs, eq, rhs_bits);
  if (!mask) return rhs == 0;
  int pivot = 63 - __clzll(mask);
  if (trail_size >= 48) return false;
  eq[pivot] = mask;
  if (rhs) rhs_bits |= 1ULL << pivot;
  else rhs_bits &= ~(1ULL << pivot);
  trail[trail_size++] = uint8_t(pivot);
  return true;
}

__device__ __forceinline__ void restore(int mark, uint64_t *eq,
                                         uint64_t &rhs_bits,
                                         uint8_t *trail,
                                         int &trail_size) {
  while (trail_size > mark) {
    int pivot = trail[--trail_size];
    eq[pivot] = 0;
    rhs_bits &= ~(1ULL << pivot);
  }
}

__global__ void solve_all(const Pair *pairs, int count,
                          const RankTwo *rank_two, int rank_two_count,
                          uint8_t *status, uint64_t *node_count,
                          uint32_t *constraint_count) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= count) return;
  const unsigned l1 = pairs[tid].a, l2 = pairs[tid].b;

  uint64_t clauses_a[192], clauses_b[192];
  int nclauses = 0;
  for (int i = 0; i < rank_two_count; ++i) {
    unsigned m = rank_two[i].message;
    if ((__popc(m & l1) & 1) || (__popc(m & l2) & 1)) continue;
    if (nclauses >= 192) { status[tid] = 2; return; }
    uint64_t a = 0, b = 0;
    for (int component = 0; component < 4; ++component) {
      if ((rank_two[i].n1 >> component) & 1) a |= uint64_t(m) << (12 * component);
      if ((rank_two[i].n2 >> component) & 1) b |= uint64_t(m) << (12 * component);
    }
    if (a > b) { uint64_t t = a; a = b; b = t; }
    clauses_a[nclauses] = a;
    clauses_b[nclauses] = b;
    ++nclauses;
  }
  constraint_count[tid] = nclauses;

  uint64_t eq[48] = {};
  uint64_t rhs_bits = 0;
  uint8_t trail[48];
  int trail_size = 0;

  uint64_t frame_a[49], frame_b[49];
  uint8_t frame_va[49], frame_vb[49], frame_mark[49], frame_stage[49];
  int depth = 0;
  uint64_t nodes = 0;

  while (true) {
    ++nodes;
    bool conflict = false, have = false;
    uint64_t chosen_a = 0, chosen_b = 0;
    int chosen_va = 0, chosen_vb = 0, best_n = 99, best_w = 999;
    for (int i = 0; i < nclauses; ++i) {
      uint64_t a = clauses_a[i], b = clauses_b[i];
      int va = 0, vb = 0;
      reduce_form(a, va, eq, rhs_bits);
      reduce_form(b, vb, eq, rhs_bits);
      if ((!a && va) || (!b && vb)) continue;
      if (!a && !va && !b && !vb) { conflict = true; break; }
      int unresolved = int(a != 0) + int(b != 0);
      int weight = __popcll(a) + __popcll(b);
      if (!have || unresolved < best_n || (unresolved == best_n && weight < best_w)) {
        have = true; best_n = unresolved; best_w = weight;
        chosen_a = a; chosen_b = b; chosen_va = va; chosen_vb = vb;
      }
    }
    if (!conflict && !have) {
      status[tid] = 1; node_count[tid] = nodes; return;
    }
    if (!conflict) {
      if (depth >= 49) { status[tid] = 3; node_count[tid] = nodes; return; }
      frame_a[depth] = chosen_a; frame_b[depth] = chosen_b;
      frame_va[depth] = uint8_t(chosen_va); frame_vb[depth] = uint8_t(chosen_vb);
      frame_mark[depth] = uint8_t(trail_size);
      if (!chosen_a) {
        frame_stage[depth++] = 2;
        if (add_equation(chosen_b, 1 ^ chosen_vb, eq, rhs_bits, trail, trail_size)) continue;
      } else if (!chosen_b) {
        frame_stage[depth++] = 2;
        if (add_equation(chosen_a, 1 ^ chosen_va, eq, rhs_bits, trail, trail_size)) continue;
      } else {
        frame_stage[depth++] = 1;
        if (add_equation(chosen_a, 1 ^ chosen_va, eq, rhs_bits, trail, trail_size)) continue;
      }
    }

    bool descended = false;
    while (depth > 0) {
      int f = depth - 1;
      restore(frame_mark[f], eq, rhs_bits, trail, trail_size);
      if (frame_stage[f] == 1) {
        frame_stage[f] = 2;
        bool ok = add_equation(frame_a[f], frame_va[f], eq, rhs_bits, trail, trail_size);
        if (ok) ok = add_equation(frame_b[f], 1 ^ frame_vb[f], eq, rhs_bits, trail, trail_size);
        if (ok) { descended = true; break; }
        restore(frame_mark[f], eq, rhs_bits, trail, trail_size);
      }
      --depth;
    }
    if (descended) continue;
    status[tid] = 0; node_count[tid] = nodes; return;
  }
}

static int rank4(uint16_t m) {
  std::array<unsigned,4> rows{};
  for (int i=0;i<4;++i) rows[i]=(m>>(4*i))&15;
  int rank=0;
  for(int bit=0;bit<4;++bit){int p=-1;for(int i=rank;i<4;++i)if((rows[i]>>bit)&1){p=i;break;}if(p<0)continue;std::swap(rows[rank],rows[p]);for(int i=rank+1;i<4;++i)if((rows[i]>>bit)&1)rows[i]^=rows[rank];++rank;}
  return rank;
}
static std::vector<unsigned> vector_basis(const std::vector<unsigned>& vectors) {
  std::array<unsigned,12> pivots{};
  for(unsigned v:vectors)while(v){int p=31-std::countl_zero(v);if(pivots[p])v^=pivots[p];else{pivots[p]=v;break;}}
  std::vector<unsigned> result;for(int p=11;p>=0;--p)if(pivots[p])result.push_back(pivots[p]);return result;
}
static std::array<unsigned,4> rows_of(uint16_t m){return {unsigned(m&15),unsigned((m>>4)&15),unsigned((m>>8)&15),unsigned((m>>12)&15)};}
static std::array<uint8_t,2> null_basis(uint16_t matrix){auto rows=rows_of(matrix);std::vector<unsigned>null;for(unsigned v=1;v<16;++v){bool ok=true;for(auto r:rows)if(std::popcount(r&v)&1){ok=false;break;}if(ok)null.push_back(v);}auto b=vector_basis(null);if(b.size()!=2){std::cerr<<"bad null\n";std::exit(2);}return{uint8_t(b[0]),uint8_t(b[1])};}

static constexpr std::array<uint16_t,12> D0={33825,17048,10692,40034,35397,42200,19754,53988,36015,51832,42828,29802};
static constexpr std::array<uint16_t,12> D2={32771,16398,8196,4104,2058,1029,524,260,140,72,47,22};
static constexpr std::array<uint16_t,12> D3={32781,16396,8206,4106,2054,1029,516,268,140,70,41,20};
static uint16_t matrix_of(unsigned x,int cls){const auto&b=cls==0?D0:(cls==2?D2:D3);uint16_t m=0;for(int i=0;i<12;++i)if((x>>i)&1)m^=b[i];return m;}

static void cuda_check(cudaError_t error,const char*where){if(error!=cudaSuccess){std::cerr<<where<<": "<<cudaGetErrorString(error)<<"\n";std::exit(3);}}
static void fnv_byte(uint64_t &hash,uint8_t value){hash^=value;hash*=1099511628211ULL;}
template<class T> static void fnv_little(uint64_t &hash,T value){for(size_t i=0;i<sizeof(T);++i)fnv_byte(hash,uint8_t((value>>(8*i))&0xff));}

int main(int argc,char**argv){
  if(argc<5){std::cerr<<"usage: class(0|2|3) start count output.json\n";return 1;}
  int cls=std::stoi(argv[1]);uint64_t start=std::stoull(argv[2]),requested=std::stoull(argv[3]);std::string output=argv[4];if(cls!=0&&cls!=2&&cls!=3)return 1;
  std::vector<Pair> all;all.reserve(2794155);for(unsigned a=1;a<4096;++a)for(unsigned b=a+1;b<4096;++b)if((a^b)>b)all.push_back({uint16_t(a),uint16_t(b)});if(all.size()!=2794155||start>=all.size())return 2;uint64_t count=std::min<uint64_t>(requested,all.size()-start);
  std::vector<RankTwo> r2;for(unsigned x=1;x<4096;++x){uint16_t m=matrix_of(x,cls);if(rank4(m)==2){auto n=null_basis(m);r2.push_back({uint16_t(x),n[0],n[1]});}}if(r2.size()!=525)return 2;
  Pair *dp=nullptr;RankTwo*dr=nullptr;uint8_t*ds=nullptr;uint64_t*dn=nullptr;uint32_t*dc=nullptr;
  cuda_check(cudaMalloc(&dp,count*sizeof(Pair)),"malloc pairs");cuda_check(cudaMalloc(&dr,r2.size()*sizeof(RankTwo)),"malloc rank2");cuda_check(cudaMalloc(&ds,count),"malloc status");cuda_check(cudaMalloc(&dn,count*sizeof(uint64_t)),"malloc nodes");cuda_check(cudaMalloc(&dc,count*sizeof(uint32_t)),"malloc constraints");
  cuda_check(cudaMemcpy(dp,all.data()+start,count*sizeof(Pair),cudaMemcpyHostToDevice),"copy pairs");cuda_check(cudaMemcpy(dr,r2.data(),r2.size()*sizeof(RankTwo),cudaMemcpyHostToDevice),"copy rank2");
  auto before=std::chrono::steady_clock::now();solve_all<<<(count+127)/128,128>>>(dp,count,dr,r2.size(),ds,dn,dc);cuda_check(cudaGetLastError(),"launch");cuda_check(cudaDeviceSynchronize(),"sync");double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
  std::vector<uint8_t>s(count);std::vector<uint64_t>nodes(count);std::vector<uint32_t>constraints(count);cuda_check(cudaMemcpy(s.data(),ds,count,cudaMemcpyDeviceToHost),"copy status");cuda_check(cudaMemcpy(nodes.data(),dn,count*sizeof(uint64_t),cudaMemcpyDeviceToHost),"copy nodes");cuda_check(cudaMemcpy(constraints.data(),dc,count*sizeof(uint32_t),cudaMemcpyDeviceToHost),"copy constraints");
  std::array<uint64_t,4>status_counts{};uint64_t sum_nodes=0,max_nodes=0,sum_constraints=0,digest=14695981039346656037ULL;uint32_t minc=999,maxc=0;for(size_t i=0;i<count;++i){if(s[i]<4)++status_counts[s[i]];sum_nodes+=nodes[i];sum_constraints+=constraints[i];max_nodes=std::max(max_nodes,nodes[i]);minc=std::min(minc,constraints[i]);maxc=std::max(maxc,constraints[i]);fnv_byte(digest,s[i]);fnv_little(digest,constraints[i]);fnv_little(digest,nodes[i]);}
  cudaDeviceProp prop{};cudaGetDeviceProperties(&prop,0);char digest_hex[17];std::snprintf(digest_hex,sizeof(digest_hex),"%016llx",(unsigned long long)digest);std::ofstream f(output);f<<"{\n  \"class\": "<<cls<<",\n  \"start\": "<<start<<",\n  \"count\": "<<count<<",\n  \"first_span\": ["<<all[start].a<<", "<<all[start].b<<", "<<(all[start].a^all[start].b)<<"],\n  \"last_span\": ["<<all[start+count-1].a<<", "<<all[start+count-1].b<<", "<<(all[start+count-1].a^all[start+count-1].b)<<"],\n  \"status_counts\": {\"UNSAT\": "<<status_counts[0]<<", \"SAT\": "<<status_counts[1]<<", \"constraint_overflow\": "<<status_counts[2]<<", \"solver_overflow\": "<<status_counts[3]<<"},\n  \"sum_nodes\": "<<sum_nodes<<",\n  \"max_nodes\": "<<max_nodes<<",\n  \"sum_constraints\": "<<sum_constraints<<",\n  \"constraint_count_range\": ["<<minc<<", "<<maxc<<"],\n  \"result_digest_fnv1a64\": \""<<digest_hex<<"\",\n  \"digest_record_format\": \"status:u8 || constraint_count:u32le || node_count:u64le, in kernel order\",\n  \"elapsed_seconds\": "<<seconds<<",\n  \"device\": \""<<prop.name<<"\"\n}\n";
  std::cout<<"class "<<cls<<" start "<<start<<" count "<<count<<" UNSAT "<<status_counts[0]<<" SAT "<<status_counts[1]<<" errors "<<(status_counts[2]+status_counts[3])<<" nodes "<<sum_nodes<<" max "<<max_nodes<<" digest "<<digest_hex<<" seconds "<<seconds<<"\n";
  cudaFree(dp);cudaFree(dr);cudaFree(ds);cudaFree(dn);cudaFree(dc);
}
