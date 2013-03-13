__constant int Iterations = 15;
__constant float MinRad2 = 0.3492f;
__constant float Scale = 2.04348f;

__constant float3 RotVector = (float3)(1.0f, 1.0f, 0.0f);
__constant float RotAngle = 0.0f;

__constant float3 Trans = (float3)(0.0365f, -1.8613f, 0.0365f);
__constant float3 Julia = (float3)(-0.6691f, -1.3028f, -0.45775f);

__constant float Epsilon = 1e-3;

float16 rotation_matrix(const float3 v, const float ang) {
  float c = cos(ang);
  float mc = 1.0f - c;
  float s = sin(ang);

  return (float16)(
    mc * v.x * v.x + c,       mc * v.x * v.y - s * v.z, mc * v.x * v.z + s * v.y, 0.0f,
    mc * v.x * v.y + s * v.z, mc * v.y * v.y + c,       mc * v.y * v.z - s * v.x, 0.0f,
    mc * v.x * v.z - s * v.y, mc * v.y * v.z + s * v.x, mc * v.z * v.z + c,       0.0f,
    0.0f,                     0.0f,                     0.0f,                     1.0f
  );
}

float3 transform(const float16 *mat, const float3 v) {
  return (float3)(dot(mat->s012, v), dot(mat->s456, v), dot(mat->s89A, v));
}

float DE(float3 pos) {
  float16 rot = rotation_matrix(normalize(RotVector), RotAngle);
  float4 scale = (float4)(Scale, Scale, Scale, fabs(Scale)) / MinRad2;
  float absScalem1 = fabs(Scale - 1.0f);
  float AbsScaleRaisedTo1mIters = pow(fabs(Scale), (float)(1-Iterations));
  float4 p = (float4)(pos, 1.0f);
  float4 p0 = (float4)(Julia, 1.0f);

  for (int i = 0; i < Iterations; i++) {
    p.xyz = transform(&rot, p.xyz);
    p.xyz = fabs(p.xyz) + Trans;
    float r2 = dot(p.xyz, p.xyz);
    p *= clamp(max(MinRad2 / r2, MinRad2), 0.0f, 1.0f);
    p = p*scale + p0;
  }
  return ((length(p.xyz) - absScalem1) / p.w - AbsScaleRaisedTo1mIters);
}

bool inside(float3 pos) {
  return DE(pos) < Epsilon;
}

#include <core.cl>
