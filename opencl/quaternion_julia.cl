#if !defined(ITERATIONS)
#  define ITERATIONS 5
#endif

__constant int Iterations = ITERATIONS;
__constant float4 C = (float4)(0.18f, 0.88f, 0.24f, 0.16f);
__constant float Threshold = 10.0f;

bool inside(float3 pos) {
  float4 p = (float4)(pos, 0.0);
  for (int i = 0; i < Iterations; i++) {
    p = (float4)(p.x*p.x-dot(p.yzw, p.yzw), (float3)(2.0f*p.x*p.yzw)) + C;
    float p2 = dot(p,p);
    if (p2 > Threshold) return false;
  }
  return true;
}

#include <core.cl>
