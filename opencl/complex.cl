float2 cMul(float2 a, float2 b);
float2 cPower(float2 z, float n);
float2 cInverse(float2 a);
float2 cDiv(float2 a, float2 b);
float2 cExp(float2 z);
float2 cLog(float2 a);
float2 cSqr(float2 z);
float2 cSin(float2 z);
float2 cCos(float2 z);
float2 cPower2(float2 z, float2 a);

float2 cMul(float2 a, float2 b) {
  return (float2)(a.x * b.x -  a.y * b.y, a.x * b.y + a.y * b.x);
}

float2 cPower(float2 z, float n) {
  float r2 = dot(z,z);
  return pow(r2, n / 2.0f) * (float2)(cos(n * atan2(z.y, z.x)), sin(n*atan(z.y/z.x)));
}

float2 cInverse(float2 a) {
  return (float2)(a.x, -a.y) / dot(a, a);
}

float2 cDiv(float2 a, float2 b) {
  return cMul(a, cInverse(b));
}

float2 cExp(float2 z) {
  return (float2)(exp(z.x) * cos(z.y), exp(z.x) * sin(z.y));
}

float2 cLog(float2 a) {
  float b =  atan2(a.y, a.x);
  if (b > 0.0) b -= 2.0f * M_PI;
  return (float2)(log(length(a)), b);
}

float2 cSqr(float2 z) {
  return (float2)(z.x * z.x - z.y * z.y, 2.0f * z.x * z.y);
}

float2 cSin(float2 z) {
  return (float2)(sin(z.x) * cosh(z.y), cos(z.x) * sinh(z.y));
}

float2 cCos(float2 z) {
  return (float2)(cos(z.x) * cosh(z.y), -sin(z.x) * sinh(z.y));
}

float2 cPower2(float2 z, float2 a) {
  return cExp(cMul(cLog(z), a));
}
