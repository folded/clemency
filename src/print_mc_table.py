# convert from pbourke's mc table to mine.
import sys

sys.stdout.write('#include "mc_table.hpp"\n\nuint8_t mc_table[256][16] = {\n')

edgecodes = [
  0x01,
  0x13,
  0x23,
  0x02,
  0x45,
  0x57,
  0x67,
  0x46,
  0x04,
  0x15,
  0x37,
  0x26
]

def tri((a, b, c)):
  if (a, b, c) == (-1, -1, -1):
    return '0x00,0x00,0x00'
  a, b, c = a, c, b
  if b < a and b < c:
    a, b, c = b, c, a
  elif c < a and c < b:
    a, b, c = c, a, b
  return '0x%02x,0x%02x,0x%02x' % ( edgecodes[a], edgecodes[b], edgecodes[c] )

rows = open('mc_base').readlines()

for i in range(256):
  rownum = (i & 0x33) | ((i & 0x88) >> 1) | ((i & 0x44) << 1)
  row = rows[rownum]
  tri_verts = map(int, row.split())
  tri_verts = [ tri_verts[x:x+3] for x in range(0, 15, 3) ]
  n_tri = len([ x for x in tri_verts if x != [ -1, -1, -1 ] ])
  args = ([ n_tri ] + [ tri(tv) for tv in tri_verts ] + [ i ])
  sys.stdout.write('  { %d, %s, %s, %s, %s, %s }, // 0x%02x\n' % tuple(args))

sys.stdout.write('};\n');
