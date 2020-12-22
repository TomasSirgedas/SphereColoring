#include "DataTypes.h"


Matrix4x4 Matrix4x4::inverted() const
{
   Matrix4x4 ret;
   ret.m[0].x = m[2].y*m[3].z*m[1].w - m[3].y*m[2].z*m[1].w + m[3].y*m[1].z*m[2].w - m[1].y*m[3].z*m[2].w - m[2].y*m[1].z*m[3].w + m[1].y*m[2].z*m[3].w;
   ret.m[1].x = m[3].x*m[2].z*m[1].w - m[2].x*m[3].z*m[1].w - m[3].x*m[1].z*m[2].w + m[1].x*m[3].z*m[2].w + m[2].x*m[1].z*m[3].w - m[1].x*m[2].z*m[3].w;
   ret.m[2].x = m[2].x*m[3].y*m[1].w - m[3].x*m[2].y*m[1].w + m[3].x*m[1].y*m[2].w - m[1].x*m[3].y*m[2].w - m[2].x*m[1].y*m[3].w + m[1].x*m[2].y*m[3].w;
   ret.m[3].x = m[3].x*m[2].y*m[1].z - m[2].x*m[3].y*m[1].z - m[3].x*m[1].y*m[2].z + m[1].x*m[3].y*m[2].z + m[2].x*m[1].y*m[3].z - m[1].x*m[2].y*m[3].z;
   ret.m[0].y = m[3].y*m[2].z*m[0].w - m[2].y*m[3].z*m[0].w - m[3].y*m[0].z*m[2].w + m[0].y*m[3].z*m[2].w + m[2].y*m[0].z*m[3].w - m[0].y*m[2].z*m[3].w;
   ret.m[1].y = m[2].x*m[3].z*m[0].w - m[3].x*m[2].z*m[0].w + m[3].x*m[0].z*m[2].w - m[0].x*m[3].z*m[2].w - m[2].x*m[0].z*m[3].w + m[0].x*m[2].z*m[3].w;
   ret.m[2].y = m[3].x*m[2].y*m[0].w - m[2].x*m[3].y*m[0].w - m[3].x*m[0].y*m[2].w + m[0].x*m[3].y*m[2].w + m[2].x*m[0].y*m[3].w - m[0].x*m[2].y*m[3].w;
   ret.m[3].y = m[2].x*m[3].y*m[0].z - m[3].x*m[2].y*m[0].z + m[3].x*m[0].y*m[2].z - m[0].x*m[3].y*m[2].z - m[2].x*m[0].y*m[3].z + m[0].x*m[2].y*m[3].z;
   ret.m[0].z = m[1].y*m[3].z*m[0].w - m[3].y*m[1].z*m[0].w + m[3].y*m[0].z*m[1].w - m[0].y*m[3].z*m[1].w - m[1].y*m[0].z*m[3].w + m[0].y*m[1].z*m[3].w;
   ret.m[1].z = m[3].x*m[1].z*m[0].w - m[1].x*m[3].z*m[0].w - m[3].x*m[0].z*m[1].w + m[0].x*m[3].z*m[1].w + m[1].x*m[0].z*m[3].w - m[0].x*m[1].z*m[3].w;
   ret.m[2].z = m[1].x*m[3].y*m[0].w - m[3].x*m[1].y*m[0].w + m[3].x*m[0].y*m[1].w - m[0].x*m[3].y*m[1].w - m[1].x*m[0].y*m[3].w + m[0].x*m[1].y*m[3].w;
   ret.m[3].z = m[3].x*m[1].y*m[0].z - m[1].x*m[3].y*m[0].z - m[3].x*m[0].y*m[1].z + m[0].x*m[3].y*m[1].z + m[1].x*m[0].y*m[3].z - m[0].x*m[1].y*m[3].z;
   ret.m[0].w = m[2].y*m[1].z*m[0].w - m[1].y*m[2].z*m[0].w - m[2].y*m[0].z*m[1].w + m[0].y*m[2].z*m[1].w + m[1].y*m[0].z*m[2].w - m[0].y*m[1].z*m[2].w;
   ret.m[1].w = m[1].x*m[2].z*m[0].w - m[2].x*m[1].z*m[0].w + m[2].x*m[0].z*m[1].w - m[0].x*m[2].z*m[1].w - m[1].x*m[0].z*m[2].w + m[0].x*m[1].z*m[2].w;
   ret.m[2].w = m[2].x*m[1].y*m[0].w - m[1].x*m[2].y*m[0].w - m[2].x*m[0].y*m[1].w + m[0].x*m[2].y*m[1].w + m[1].x*m[0].y*m[2].w - m[0].x*m[1].y*m[2].w;
   ret.m[3].w = m[1].x*m[2].y*m[0].z - m[2].x*m[1].y*m[0].z + m[2].x*m[0].y*m[1].z - m[0].x*m[2].y*m[1].z - m[1].x*m[0].y*m[2].z + m[0].x*m[1].y*m[2].z;
   double det     = m[0].x* ret.m[0].x + m[1].x* ret[0].y  + m[2].x* ret[0].z + m[3].x* ret[0].w;
   double invdet  = 1 / det;
   ret.m[0] *= invdet;
   ret.m[1] *= invdet;
   ret.m[2] *= invdet;
   ret.m[3] *= invdet;
   return ret;
}

Matrix4x4 Matrix4x4::pow( int pwr ) const
{
   Matrix4x4 ret;
   for ( int i = 0; i < abs(pwr); i++ )
      ret = ret * *this;
   return pwr > 0 ? ret : ret.inverted();
}