#include "Simulation.h"
#include <set>
#include <map>

using namespace std;

void Simulation::init( shared_ptr<Dual> dual, shared_ptr<Graph> graph, double radius )
{
   _Dual = dual;
   _Graph = graph;

   if ( _Graph )
   {
      _KeepCloseFars = _Graph->calcKeepCloseFars();
      _LineVertexConstraints = _Graph->calcLineVertexConstraints();
   }
   
   //for ( const Graph::KeepCloseFar& kcf : _KeepCloseFars )
   //   if ( kcf.keepClose )
   //      qDebug() << graph->idOf( kcf.a ) << "-" << graph->idOf( kcf.b );

   _Radius = radius;
   normalizeVertices();
}


void Simulation::normalizeVertices()
{
   for ( Dual::Vertex& vtx : _Dual->_Vertices )
      vtx._Pos = vtx._Pos.normalized() * _Radius;

   if ( !_Graph )
      return;
   for ( Graph::Vertex& vtx : _Graph->_Vertices )
      vtx._Pos = vtx._Pos.normalized() * _Radius;
}

double Simulation::step( double& paddingError )
{   
   if ( !_Graph )
      return -1;
   static bool s_printErrors = false;
   bool printErrors = s_printErrors && _PaddingError == 0;
   double totalError = 0;
   paddingError = 0;
   vector<XYZ> vel( _Graph->_Vertices.size() );

   for ( const Graph::KeepCloseFar& kcf : _KeepCloseFars )
   {
      XYZ a = _Graph->posOf( kcf.a );
      XYZ b = _Graph->posOf( kcf.b );
      double dist = a.dist( b );

      double pad = kcf.keepClose && kcf.keepFar ? 0 : _Padding;
      if ( kcf.keepClose && dist >= 1.-pad )
      {
         vel[kcf.a._Index] += (kcf.a._Mtx.inverted() * (b-a)) * (dist-(1-pad)) * .03;
         vel[kcf.b._Index] += (kcf.b._Mtx.inverted() * (a-b)) * (dist-(1-pad)) * .03;
         totalError += max(0.,dist-1);
         paddingError += dist-(1-pad);
         if ( printErrors && !kcf.keepFar && dist-1 > 0 ) qDebug() << "keep close" << kcf.a._Index << kcf.b._Index << dist-1;
      }
      if ( kcf.keepFar && dist <= 1.+pad )
      {
         vel[kcf.a._Index] += (kcf.a._Mtx.inverted() * (a-b).normalized()) * ((1+pad)-dist) * .03;
         vel[kcf.b._Index] += (kcf.b._Mtx.inverted() * (b-a).normalized()) * ((1+pad)-dist) * .03;
         totalError += max(0.,1-dist);
         paddingError += (1+pad)-dist;
         if ( printErrors && !kcf.keepClose && 1-dist > 0 ) qDebug() << "keep far" << kcf.a._Index << kcf.b._Index << 1-dist;
      }
   }
   //static bool s_dolvc = true;
   for ( const Graph::LineVertexConstraint& lvc : _LineVertexConstraints )  if ( !lvc.curveCenter.isValid() )
   {
      double pad = _Padding;

      // if the curveCenter is valid, then the curve is a circle with radius=1, centered at lvc.curveCenter
      // otherwise,                        the curve is a circle with radius=_Radius, centered at the origin
      XYZ center = lvc.curveCenter.isValid() ? _Graph->posOf( lvc.curveCenter ) : XYZ();
      XYZ b = _Graph->posOf( lvc.b ) - center;
      XYZ a0 = _Graph->posOf( lvc.a0 ) - center;
      XYZ a1 = _Graph->posOf( lvc.a1 ) - center;
      double R = lvc.curveCenter.isValid() ? 1 : _Radius;
      double dot = a0 * a1;
      // b = x * a0 + y * a1 + z * (a0^a1)   // solve for x/y/z
      double x = (b*a0 * R*R - b*a1 * dot) / (R*R*R*R - dot*dot);
      double y = (b*a1 * R*R - b*a0 * dot) / (R*R*R*R - dot*dot);
      if ( x < 0 || x > 1 || y < 0 || y > 1 )
         continue; // b doesn't lie between a0, a1  (so the distance checks on the endpoints are sufficient)
      XYZ q = (a0*x + a1*y).normalized() * R; // project the point to the circle 
      double dist = q.dist( b );
      if ( dist >= 1+pad ) 
         continue;
      
      XYZ qb = (q-b).normalized();
      vel[lvc.a0._Index] += (lvc.a0._Mtx.inverted() * qb) * ((1+pad)-dist) *  .005;
      vel[lvc.a1._Index] += (lvc.a1._Mtx.inverted() * qb) * ((1+pad)-dist) *  .005;
      vel[lvc.b._Index]  += (lvc.b._Mtx.inverted()  * qb) * ((1+pad)-dist) * -.01;
      totalError += max(0.,1-dist);
      paddingError += (1+pad)-dist;
      if ( printErrors && 1-dist > 0 ) qDebug() << "straight line to vertex" << lvc.a0._Index << lvc.a1._Index << lvc.b._Index << 1-dist;
   }
   for ( const Graph::LineVertexConstraint& lvc : _LineVertexConstraints )  if ( lvc.curveCenter.isValid() )
   {
      double pad = _Padding;
      
      // if the curveCenter is valid, then the curve is a circle with radius=1, centered at lvc.curveCenter
      // otherwise,                        the curve is a circle with radius=_Radius, centered at the origin
      XYZ center = lvc.curveCenter.isValid() ? _Graph->posOf( lvc.curveCenter ) : XYZ();
      XYZ b = _Graph->posOf( lvc.b ) - center;
      if ( b.len2() >= (2+pad)*(2+pad) )
         continue;
      XYZ a0 = (_Graph->posOf( lvc.a0 ) - center).normalized();
      XYZ a1 = (_Graph->posOf( lvc.a1 ) - center).normalized();
      double R = lvc.curveCenter.isValid() ? 1 : _Radius;
      double dot = a0 * a1;
      // b = x * a0 + y * a1 + z * (a0^a1)   // solve for x/y/z
      double x = (b*a0 * 1 - b*a1 * dot) / (1 - dot*dot);
      double y = (b*a1 * 1 - b*a0 * dot) / (1 - dot*dot);
      if ( x < 0 || y < 0 )
         continue; // b doesn't lie between a0, a1  (so the distance checks on the endpoints are sufficient)
      XYZ q = (a0*x + a1*y).normalized() * R; // project the point to the circle 
      double dist = q.dist( b );
      if ( dist >= 1+pad ) 
         continue;
      if ( isnan(dist) )
         qDebug() << lvc.curveCenter._Index << lvc.a0._Index << lvc.a1._Index << lvc.b._Index << dist;

      XYZ qb = (q-b).normalized();
      vel[lvc.curveCenter._Index]  += (lvc.curveCenter._Mtx.inverted() * qb) * ((1+pad)-dist) *  .00003;
      vel[lvc.a0._Index]           += (lvc.a0._Mtx.inverted()          * qb) * ((1+pad)-dist) *  .00003;
      vel[lvc.a1._Index]           += (lvc.a1._Mtx.inverted()          * qb) * ((1+pad)-dist) *  .00003;
      vel[lvc.b._Index]            += (lvc.b._Mtx.inverted()           * qb) * ((1+pad)-dist) * -.00009;
      totalError += max(0.,1-dist);
      paddingError += (1+pad)-dist;
      if ( printErrors && 1-dist > 0 ) qDebug() << "curved line to vertex" << lvc.a0._Index << lvc.a1._Index << lvc.curveCenter._Index << lvc.b._Index << 1-dist;
   }


   // apply velocities
   for ( int i = 0; i < (int)vel.size(); i++ ) if ( !_Graph->_Vertices[i]._IsSymmetrical )
   {
      _Graph->_Vertices[i]._Pos += vel[i];
   }
   normalizeVertices();

   return totalError;
}

double Simulation::step( int numSteps )
{
   double tot = 0;
   double totalPaddingError = 0;
   for ( int i = 0; i < numSteps; i++ )
   {
      double paddingError = 0;
      tot += step( paddingError );
      totalPaddingError += paddingError;
   }

   //{
   //   XYZ a = _Graph->posOf( _Graph->fromId( 6 ) );
   //   XYZ b = _Graph->posOf( _Graph->fromId( 8 ) );
   //   qDebug() << "dist = " << a.dist( b );
   //}


   _PaddingError = totalPaddingError / numSteps;
   return tot / numSteps;
}