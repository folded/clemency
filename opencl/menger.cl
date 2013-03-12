#if !defined(ITERATIONS)
#  define ITERATIONS 7
#endif
__constant int Iterations = ITERATIONS;

bool inside(float3 pos) {
  if(pos.x < 0.0f || pos.x > 1.0f ||
     pos.y < 0.0f || pos.y > 1.0f ||
     pos.z < 0.0f || pos.z > 1.0f) return false;

  for (int m = 0; m < Iterations; m++) {
    pos *= 3.0f;

    if ((pos.x > 1.0 && pos.x < 2.0 && pos.y > 1.0 && pos.y < 2.0) ||
        (pos.y > 1.0 && pos.y < 2.0 && pos.z > 1.0 && pos.z < 2.0) ||
        (pos.x > 1.0 && pos.x < 2.0 && pos.z > 1.0 && pos.z < 2.0)) {
      return false;
    }

    pos = fmod(pos, (float3)(1.0));
  }
  return true;
}

#include <core.cl>
