#pragma once

#include <vector>
#include <QDebug>
#include <QMatrix4x4>
#include <QPolygon>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include "DataTypes.h"

using namespace std;


static int mod( int x, int m ) { return (x%m+m)%m; }

class Perm
{
public:
   Perm() {}
   Perm( int n ) { for ( int i = 0; i < n; i++ ) v.push_back( i ); }
   Perm( const vector<int>& v ) : v(v) {}
   int operator[]( int idx ) const { return v[idx]; }
   int& operator[]( int idx ) { return v[idx]; }

   Perm operator*( const Perm& rhs ) const
   {
      Perm ret;
      for ( int i = 0; i < (int) v.size(); i++ )
         ret.v.push_back( rhs[v[i]] );
      return ret;
   }
   Perm inverted() const
   {
      Perm ret( v.size() );
      for ( int i = 0; i < (int) v.size(); i++ )
         ret[v[i]] = i;
      return ret;
   }
   Perm pow( int n ) const
   {
      Perm ret( v.size() );
      for ( int i = 0; i < n; i++ )
         ret = *this * ret;
      return ret;
   }
   bool circularlyMatches( const Perm& rhs ) const
   {
      int idxOffset = rhs.inverted()[v[0]];
      for ( int i = 0; i < (int) rhs.v.size(); i++ )
         if ( v[i] != rhs[mod(i+idxOffset,(int)v.size())] )
            return false;
      return true;
   }

public:
   vector<int> v;
};
static QDebug operator<<( const QDebug& d, const Perm& perm ) { return d << perm.v; }



const double PHI = .5 + sqrt(1.25);

//typedef QVector3D XYZ;
//typedef QMatrix4x4 QMtx4x4;
typedef Matrix4x4 QMtx4x4;



bool isClockwiseTri( const QPolygonF& poly );

QMtx4x4 toMatrix( const XYZ& a, const XYZ& b, const XYZ& c );
QMtx4x4 translation( const XYZ& p );
QMtx4x4 map( const vector<XYZ>& a, const vector<XYZ>& b );
QMtx4x4 matrixRotateToZAxis( XYZ& p );
QMtx4x4 pow( const QMtx4x4& m, int power );
bool fuzzyCompare( const QMtx4x4& a, const QMtx4x4& b );
bool isIdentity( const QMtx4x4& a );
static XYZ operator*( const QMtx4x4& m, const XYZ& p ) 
{ 
   //QVector3D ret = m * QVector3D( p.x, p.y, p.z );
   //return XYZ( ret.x(), ret.y(), ret.z() );
   return (m * XYZW( p )).toXYZ();
}
static vector<XYZ> operator*( const QMtx4x4& m, const vector<XYZ>& v ) 
{ 
   vector<XYZ> ret;
   for ( const XYZ& p : v )
      ret.push_back( m * p );
   return ret;
}

class ISymmetry
{
public:
   struct Config
   {
      QMtx4x4 m;
      vector<int> state;
      bool isHomeState() const 
      { 
         for ( int x : state )
            if ( x != 0 )
               return false;
         return true;
      }
   };

public:
   virtual Perm colorPermOf( const QMtx4x4& m ) const = 0;
   virtual vector<Config> matrices() const = 0;
};

class IcoSymmetry : public ISymmetry
{
public:
   IcoSymmetry() 
   {
      _Pts = { XYZ(    -1,    0, -PHI )  
             , XYZ(     1,    0, -PHI )  
             , XYZ(     0,  PHI,   -1 )  
             , XYZ(  -PHI,    1,    0 )  
             , XYZ(  -PHI,   -1,    0 )  
             , XYZ(     0, -PHI,   -1 )
             , XYZ(     1,    0,  PHI )  
             , XYZ(    -1,    0,  PHI )  
             , XYZ(     0, -PHI,    1 )  
             , XYZ(   PHI,   -1,    0 )  
             , XYZ(   PHI,    1,    0 )  
             , XYZ(     0,  PHI,    1 ) };
   }
   XYZ operator[]( int idx ) const { return _Pts[idx]; }
   QMtx4x4 map( const vector<int>& a, const vector<int>& b ) const
   {
      return ::map( vector<XYZ> { _Pts[a[0]], _Pts[a[1]], _Pts[a[2]] }, vector<XYZ> { _Pts[b[0]], _Pts[b[1]], _Pts[b[2]] } );
   }
   vector<Config> matrices() const
   {
      vector<Config> ret;

      for ( int sym0 : { 0, 1, 2 } )
      for ( int sym1 : { 0, 1 } )
      for ( int sym2 : { 0, 1 } )
      for ( int sym3 : { 0, 1, 2, 3, 4 } )
      {
         QMtx4x4 m = pow( map( {0,1,2}, {3,0,2} ), sym3 )
                   * pow( map( {0,1,2}, {7,6,8} ), sym2 )
                   * pow( map( {0,1,2}, {1,0,5} ), sym1 )
                   * pow( map( {0,1,2}, {1,2,0} ), sym0 );
         ret.push_back( Config { m, {sym0,sym1,sym2,sym3} } );
      }
      return ret;
   }
   int id( const XYZ& p ) const
   {
      for ( int i = 0; i < (int) _Pts.size(); i++ )
         if ( _Pts[i].dist( p ) < 1e-5 )
            return i;
      throw 777;
   }
   double radius() const { return _Pts[0].len(); }

   int tetrColorOf( const XYZ& p ) const
   {
      if ( p.dist2( _Pts[0]+_Pts[1]+_Pts[2] ) < 1e-9 ) return 0;
      if ( p.dist2( _Pts[3]+_Pts[7]+_Pts[11] ) < 1e-9 ) return 0;
      if ( p.dist2( _Pts[6]+_Pts[9]+_Pts[10] ) < 1e-9 ) return 0;
      if ( p.dist2( _Pts[4]+_Pts[5]+_Pts[8] ) < 1e-9 ) return 0;
      
      if ( p.dist2( _Pts[0]+_Pts[2]+_Pts[3] ) < 1e-9 ) return 1;
      if ( p.dist2( _Pts[6]+_Pts[10]+_Pts[11] ) < 1e-9 ) return 1;
      if ( p.dist2( _Pts[1]+_Pts[5]+_Pts[9] ) < 1e-9 ) return 1;
      if ( p.dist2( _Pts[4]+_Pts[7]+_Pts[8] ) < 1e-9 ) return 1;

      if ( p.dist2( _Pts[0]+_Pts[3]+_Pts[4] ) < 1e-9 ) return 2;
      if ( p.dist2( _Pts[1]+_Pts[2]+_Pts[10] ) < 1e-9 ) return 2;
      if ( p.dist2( _Pts[5]+_Pts[8]+_Pts[9] ) < 1e-9 ) return 2;
      if ( p.dist2( _Pts[6]+_Pts[7]+_Pts[11] ) < 1e-9 ) return 2;

      if ( p.dist2( _Pts[0]+_Pts[4]+_Pts[5] ) < 1e-9 ) return 3;
      if ( p.dist2( _Pts[1]+_Pts[9]+_Pts[10] ) < 1e-9 ) return 3;
      if ( p.dist2( _Pts[2]+_Pts[3]+_Pts[11] ) < 1e-9 ) return 3;
      if ( p.dist2( _Pts[6]+_Pts[7]+_Pts[8] ) < 1e-9 ) return 3;

      return 4;
   }

   Perm tetrColorPermOf( const QMtx4x4& m ) const
   {
      XYZ p012 = m * (_Pts[0]+_Pts[1]+_Pts[2]);
      XYZ p023 = m * (_Pts[0]+_Pts[2]+_Pts[3]);
      XYZ p034 = m * (_Pts[0]+_Pts[3]+_Pts[4]);
      XYZ p045 = m * (_Pts[0]+_Pts[4]+_Pts[5]);
      XYZ p051 = m * (_Pts[0]+_Pts[5]+_Pts[1]);

      return Perm( {tetrColorOf(p012),tetrColorOf(p023),tetrColorOf(p034),tetrColorOf(p045),tetrColorOf(p051),5,6} );  
   }

   Perm colorPermOf( const QMtx4x4& m ) const
   {            
      return Perm( { id( m*_Pts[0] )%6, id( m*_Pts[1] )%6, id( m*_Pts[2] )%6, id( m*_Pts[3] )%6, id( m*_Pts[4] )%6, id( m*_Pts[5] )%6, 6} );        
      //return tetrColorPermOf( m );
   }

public:
   vector<XYZ> _Pts;
};

class GlobalSymmetry
{
private:
   GlobalSymmetry() { _Symmetry.reset( new IcoSymmetry ); }

public:
   static GlobalSymmetry& theInstance() { static GlobalSymmetry s_theInstance; return s_theInstance; }
   static ISymmetry* symmetry() { return theInstance()._Symmetry.get(); }
   static Perm colorPermOf( const QMtx4x4& m ) { return symmetry()->colorPermOf( m ); }
   static vector<ISymmetry::Config> matrices() { return symmetry()->matrices(); }
      
public:
   shared_ptr<ISymmetry> _Symmetry;
};

uint64_t matrixId( const QMtx4x4& m );

class MatrixIndexMap
{
private:
   MatrixIndexMap()
   {
      _Matrices = GlobalSymmetry::matrices();
      for ( int i = 0; i < (int) _Matrices.size(); i++ )
         _MatrixIdToIndex[matrixId( _Matrices[i].m )] = i;
   }
public:
   static MatrixIndexMap& theInstance() { static MatrixIndexMap s_theInstance; return s_theInstance; }
   static int indexOf( const QMtx4x4& m ) 
   { 
      int id = matrixId( m );
      //assert( theInstance()._MatrixIdToIndex.count( id );
      return theInstance()._MatrixIdToIndex.at( id ); 
   }
   static QMtx4x4 at( int index ) { return theInstance()._Matrices[index].m; }

public:
   vector<ISymmetry::Config> _Matrices;
   unordered_map<uint64_t, int> _MatrixIdToIndex;
};

class MatrixSymmetryMap
{
public:
   MatrixSymmetryMap( const XYZ& symmetricalPt )
   {
      const vector<ISymmetry::Config>& m = MatrixIndexMap::theInstance()._Matrices;
      for ( int a = 0; a < (int) m.size(); a++ )
      {
         _MapToReal.push_back( a );
         for ( int i = 0; i < a; i++ )
            if ( (m[a].m * symmetricalPt).dist2( m[i].m * symmetricalPt ) < 1e-8 )
            {
               _MapToReal.back() = i;
               break;
            }
      }
      _SymmetricMatrices.resize( m.size() );
      for ( int a = 0; a < (int) m.size(); a++ )
         for ( int b = 0; b < (int) m.size(); b++ )
            if ( _MapToReal[MatrixIndexMap::indexOf( m[a].m )] == _MapToReal[MatrixIndexMap::indexOf( m[b].m )] )
               _SymmetricMatrices[a].push_back( b );
   }
   bool isReal( const QMtx4x4& m ) const { int index = MatrixIndexMap::indexOf( m ); return _MapToReal[index] == index; }
   QMtx4x4 toReal( const QMtx4x4& m ) const { return MatrixIndexMap::at( _MapToReal[MatrixIndexMap::indexOf( m )] ); }
   bool match( const QMtx4x4& a, const QMtx4x4& b ) const { return _MapToReal[MatrixIndexMap::indexOf( a )] == _MapToReal[MatrixIndexMap::indexOf( b )]; }
   vector<QMtx4x4> symmetricMatrices( const QMtx4x4& m ) const 
   { 
      vector<QMtx4x4> ret;
      for ( int idx : _SymmetricMatrices[MatrixIndexMap::indexOf( m )] ) 
         ret.push_back( MatrixIndexMap::at( idx ) ); 
      return ret;
   }

   //static MatrixSymmetryMap* symmetryNone()  { static IcoSymmetry ico; static MatrixSymmetryMap s_map( XYZ(7,8,9) ); return &s_map; }
   //static MatrixSymmetryMap* symmetry012()   { static IcoSymmetry ico; static MatrixSymmetryMap s_map( ico[0]+ico[1]+ico[2] ); return &s_map; }
   //static MatrixSymmetryMap* symmetry12345() { static IcoSymmetry ico; static MatrixSymmetryMap s_map( ico[0] ); return &s_map; }
   //static MatrixSymmetryMap* symmetry01()    { static IcoSymmetry ico; static MatrixSymmetryMap s_map( ico[0]+ico[1] ); return &s_map; }
   //static MatrixSymmetryMap* symmetryFor( const XYZ& p )
   //{
   //   IcoSymmetry ico;
   //   if ( p.normalized().dist2( ico[0].normalized() ) < 1e-9 )
   //      return symmetry12345();
   //   if ( p.normalized().dist2( ( ico[0] + ico[1] + ico[2] ).normalized() ) < 1e-9 )
   //      return symmetry012();
   //   if ( p.normalized().dist2( ( ico[0] + ico[1] ).normalized() ) < 1e-9 )
   //      return symmetry01();
   //   return symmetryNone();
   //}

   static shared_ptr<MatrixSymmetryMap> symmetryNone() 
   { 
      static shared_ptr<MatrixSymmetryMap> s_noSymmetry( new MatrixSymmetryMap(XYZ(7,8,9)) ); 
      return s_noSymmetry;
   }


   static shared_ptr<MatrixSymmetryMap> symmetryFor( const XYZ& p )
   {
      shared_ptr<MatrixSymmetryMap> ret( new MatrixSymmetryMap( p ) );
      if ( !ret->hasSymmetry() )
         return symmetryNone(); // less duplication
      return ret;
   }


   bool hasSymmetry() const { return _SymmetricMatrices[0].size() > 1; }

public:
   vector<int> _MapToReal;
   vector<vector<int>> _SymmetricMatrices; // (inverse of _MapToReal)
};

class Graph
{
public:
   class VertexPtr
   {
   public:
      VertexPtr() : _Index(-1) {}
      VertexPtr( int idx, const QMtx4x4& mtx ) : _Index(idx), _Mtx(mtx) {}
      bool isValid() const { return _Index >= 0; }
      VertexPtr premul( const QMtx4x4& mtx ) const { return VertexPtr( _Index, mtx * _Mtx ); }
      bool operator==( const VertexPtr& rhs ) const;
      uint64_t id() const;

      int _Index;
      QMtx4x4 _Mtx;
   };
   class TilePtr
   {
   public:
      TilePtr() : _Index(-1) {}
      TilePtr( int idx, const QMtx4x4& mtx ) : _Index(idx), _Mtx(mtx) {}
      TilePtr premul( const QMtx4x4& mtx ) const { return TilePtr( _Index, mtx * _Mtx ); }
      bool isValid() const { return _Index >= 0; }
      //bool operator==( const TilePtr& rhs ) const;
      int _Index;
      QMtx4x4 _Mtx;
   };
   class Vertex
   {
   public:
      Vertex( int idx ) : _Index(idx) {}
      int _Index;
      bool _IsSymmetrical;
      XYZ _Pos;
      vector<VertexPtr> _Neighbors;
      vector<TilePtr> _Tiles;
   };
   class Tile
   {
   public:
      int _Index;
      int _Color;
      vector<VertexPtr> _Vertices;
      //bool _IsSymmetrical = false;
      shared_ptr<MatrixSymmetryMap> _SymmetryMap;

      bool hasVertex( const VertexPtr& a ) const { 
         for ( const VertexPtr& b : _Vertices )
            if ( a._Index == b._Index && fuzzyCompare( a._Mtx, b._Mtx ) )
               return true;
         return false;
      }
   };
   struct KeepCloseFar
   {
      VertexPtr a;
      VertexPtr b;
      bool keepClose;
      bool keepFar;
   };
   // these two shouldn't get too close together
   // - line[a0,a1] curves centered on curveCenter
   // - vertex b
   struct LineVertexConstraint
   {
      VertexPtr a0;  
      VertexPtr a1;
      VertexPtr curveCenter;
      VertexPtr b;
   };

public:
   Graph();
   XYZ posOf( const VertexPtr& vtx ) const;
   //XYZ originalPosOf( const VertexPtr& vtx ) const;
   int idOf( const VertexPtr& vtx ) const;
   VertexPtr fromId( int id ) const;
   void addNeighbor( const VertexPtr& a, const VertexPtr& b );
   vector<VertexPtr> neighbors( const VertexPtr& vtx ) const;
   vector<VertexPtr> neighbors( const VertexPtr& vtx, int depth ) const;
   VertexPtr operator[]( int idx ) const { return VertexPtr( idx, QMtx4x4() ); }
   vector<int> colorsAt( const VertexPtr& vtx ) const;
   uint32_t colorBits( const VertexPtr& vtx ) const;
   int colorOf( const TilePtr& tile ) const;
   vector<TilePtr> tilesAt( const VertexPtr& vtx ) const;
   vector<TilePtr> tilesAt( const VertexPtr& a, const VertexPtr& b ) const;
   TilePtr tileWithColor( const VertexPtr& vtx, int color ) const;
   bool mustBeFar( const VertexPtr& a, const VertexPtr& b ) const;
   bool mustBeClose( const VertexPtr& a, const VertexPtr& b ) const;
   bool eq( const TilePtr& a, const TilePtr& b ) const;
   bool eq( const VertexPtr& a, const VertexPtr& b ) const;
   vector<VertexPtr> allVertices() const;
   vector<VertexPtr> rawVertices() const;
   vector<TilePtr> rawTiles() const;
   vector<TilePtr> allTiles() const;
   vector<KeepCloseFar> calcKeepCloseFars() const;
   vector<LineVertexConstraint> calcLineVertexConstraints() const;
   vector<pair<VertexPtr,VertexPtr>> calcPerimeter() const;
   VertexPtr calcCurve( const VertexPtr& a, const VertexPtr& b ) const;
   vector<VertexPtr> verticesForTile( const TilePtr& tile ) const;
   

private:
   void neighbors( const VertexPtr& vtx, int depth, vector<VertexPtr>& v, unordered_set<uint64_t>& st ) const;

public:
   vector<Vertex> _Vertices;
   vector<Tile> _Tiles;
};


class Dual
{
public:
   class VertexPtr
   {
   public:
      VertexPtr() : _Index(-1) {}
      VertexPtr( int idx, const QMtx4x4& mtx ) : _Index(idx), _Mtx(mtx) {}
      bool isValid() const { return _Index >= 0; }
      //VertexPtr premul( const QMtx4x4& mtx ) const { return VertexPtr( _Index, mtx * _Mtx ); }
      bool operator==( const VertexPtr& rhs ) const;

      int _Index;
      QMtx4x4 _Mtx;
   };
   class Vertex
   {
   public:
      Vertex( int idx ) : _Index(idx) {}
      VertexPtr toVertexPtr() const { return VertexPtr( _Index, QMtx4x4() ); }
      int findNeighborIdx( const VertexPtr& a ) const;
      bool hasNeighbor( const VertexPtr& a ) const;
      void eraseNeighbor( const VertexPtr& a );
      int _Index;
      //bool _IsSymmetrical;
      shared_ptr<MatrixSymmetryMap> _SymmetryMap = MatrixSymmetryMap::symmetryNone();
      XYZ _Pos;
      int _Color;
      vector<VertexPtr> _Neighbors;
   };

public:
   void addVertex( int color, const XYZ& pos )
   {
      Vertex vtx( (int) _Vertices.size() );
      //vtx._IsSymmetrical = false;
      vtx._Pos = pos;
      vtx._Color = color;
      _Vertices.push_back( vtx );
   }
   vector<VertexPtr> allVertices() const;
   vector<VertexPtr> baseVertices() const;
   XYZ posOf( const VertexPtr& vtx ) const;
   void setPos( const VertexPtr& vtx, const XYZ& pos );
   int colorOf( const VertexPtr& vtx ) const;
   void setColorOf( const VertexPtr& vtx, int color );
   void toggleEdge( const VertexPtr& a, const VertexPtr& b, bool onlyAdd=false );
   vector<VertexPtr> neighborsOf( const VertexPtr& a ) const;
   vector<VertexPtr> sortedNeighborsOf( const VertexPtr& a ) const;
   bool isDuplicate( const VertexPtr& a ) const;
   VertexPtr toReal( const VertexPtr& a ) const;
   int idOf( const VertexPtr& a ) const;
   VertexPtr fromId( int id ) const;
   VertexPtr next( const VertexPtr& a, const VertexPtr& b ) const;
   vector<Dual::VertexPtr> polygon( const VertexPtr& a, const VertexPtr& b ) const;
   VertexPtr premul( const VertexPtr& vtx, const QMtx4x4& mtx ) const;
   
public:
   vector<Vertex> _Vertices;
};