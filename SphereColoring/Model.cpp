#include "Model.h"

#include <algorithm>
#include <map>


bool isClockwiseTri( const QPolygonF& poly )
{
   QPointF u = poly.at(1) - poly.at(0);
   QPointF v = poly.at(2) - poly.at(0);
   return u.x() * v.y() < u.y() * v.x();
}

QMtx4x4 toMatrix( const XYZ& a, const XYZ& b, const XYZ& c )
{
   return QMtx4x4( a.x, b.x, c.x, 0,
                   a.y, b.y, c.y, 0, 
                   a.z, b.z, c.z, 0, 
                   0,   0,   0  , 1 );
}

vector<int> canonicalRotation( vector<int> v )
{
   rotate( v.begin(), min_element( v.begin(), v.end() ), v.end() );
   return v;
}

vector<int> sorted( vector<int> v )
{
   sort( v.begin(), v.end() );
   return v;
}

uint32_t toBitFlag( const vector<int>& v )
{
   uint32_t ret = 0;
   for ( int x : v )
      ret |= 1u << x;
   return ret;
}

QMtx4x4 translation( const XYZ& p ) { QMtx4x4 m; m.translate( p.x, p.y, p.z ); return m; }

QMtx4x4 map( const vector<XYZ>& a, const vector<XYZ>& b )
{
   if ( a.size() == 3 && b.size() == 3 )
      return toMatrix( b[0], b[1], b[2] ) * toMatrix( a[0], a[1], a[2] ).inverted();
   if ( a.size() == 4 && b.size() == 4 )
   {
      QMtx4x4 m = ::map( { a[1]-a[0], a[2]-a[0], a[3]-a[0] }, { b[1]-b[0], b[2]-b[0], b[3]-b[0] } );
      return translation( b[0] ) * m * translation( -a[0] );
   }
   throw 777;
}

// rotation matrix that rotates p to the Z-axis
QMtx4x4 matrixRotateToZAxis( XYZ& p )
{
   XYZ newZ = p.normalized();
   XYZ q = abs(newZ.z) < abs(newZ.y) ? XYZ(0,0,1) : XYZ(0,1,0);
   XYZ newX = (newZ ^ q).normalized();
   XYZ newY = (newZ ^ newX).normalized();
   return QMtx4x4( toMatrix( newX, newY, newZ ) ).inverted();
}

QMtx4x4 pow( const QMtx4x4& m, int power )
{
   QMtx4x4 ret;
   for ( int i = 0; i < power; i++ )
      ret = ret * m;
   return ret;
}

bool fuzzyCompare( const QMtx4x4& a, const QMtx4x4& b )
{
   for ( int y = 0; y < 4; y++ )
      for ( int x = 0; x < 4; x++ )
         if ( abs( a(y,x) - b(y,x) ) >= 1e-5 )
            return false;
   return true;
}

bool isIdentity( const QMtx4x4& a )
{
   return fuzzyCompare( a, QMtx4x4() );
}

bool contains( const vector<QPoint>& v, const QPoint& p )
{
   return find( v.begin(), v.end(), p ) != v.end();
}

// [0..2^18)
// [0..262144)
uint64_t matrixId( const QMtx4x4& m )
{
   uint64_t ret = 0;
   ret = ret * 8 + ( int( m(0,0) / .24 ) + 3 );
   ret = ret * 8 + ( int( m(1,0) / .24 ) + 3 );
   ret = ret * 8 + ( int( m(2,0) / .24 ) + 3 );
   ret = ret * 8 + ( int( m(0,1) / .24 ) + 3 );
   ret = ret * 8 + ( int( m(1,1) / .24 ) + 3 );
   ret = ret * 8 + ( int( m(2,1) / .24 ) + 3 );
   return ret;
}




namespace
{
   bool isMissingEdge( const QPoint& a, const QPoint& b )
   {
      static vector<pair<QPoint, QPoint>> s_missingEdges = 
      { { {3,-3}, {3,-4} },
      { {3,-3}, {2,-2} },
      { {3,-2}, {4,-3} },
      { {2,-1}, {1,-1} },
      { {3,0}, {4,0} },
      { {4,-1}, {5,-1} },
      { {6,-1}, {6,0} },
      { {2,-2}, {1,-2} },
      { {2,-2}, {3,-3} } };
      for ( const auto& edge : s_missingEdges ) 
      {
         if ( edge.first == a && edge.second == b )
            return true;
         if ( edge.first == b && edge.second == a )
            return true;
      }
      return false;
   }
}


Graph::Graph()
{
}

XYZ avg( const vector<XYZ>& pts )
{
   XYZ sum(0,0,0);
   for ( const XYZ& pt : pts )
      sum += pt;
   return sum/pts.size();
}

XYZ Graph::posOf( const VertexPtr& vtx ) const
{
   return vtx._Mtx * _Vertices[vtx._Index]._Pos;
}

int Graph::idOf( const VertexPtr& vtx ) const
{
   return MatrixIndexMap::indexOf( vtx._Mtx ) * _Vertices.size() + vtx._Index;
}

Graph::VertexPtr Graph::fromId( int id ) const
{
   int sz = (int) _Vertices.size();   
   return VertexPtr( id % sz, MatrixIndexMap::at( id / sz ) );
}

uint64_t Graph::VertexPtr::id() const
{
   return _Index * (1LL<<18) + matrixId( _Mtx );
}

void Graph::addNeighbor( const VertexPtr& a, const VertexPtr& b )
{
   VertexPtr bb = b.premul( a._Mtx.inverted() );

   if ( a == b )
      qDebug( "addNeighbor dup" );
   
   for ( const VertexPtr& v : _Vertices[a._Index]._Neighbors )
      if ( bb == v )
         return; // already have it

   _Vertices[a._Index]._Neighbors.push_back( bb );
}

vector<Graph::VertexPtr> Graph::neighbors( const VertexPtr& vtx ) const
{
   vector<Graph::VertexPtr> ret;
   
   for ( const Graph::VertexPtr& neighb : _Vertices[vtx._Index]._Neighbors )
      ret.push_back( neighb.premul( vtx._Mtx ) );

   return ret;
}

vector<Graph::VertexPtr> Graph::neighbors( const VertexPtr& vtx, int depth ) const
{
   unordered_set<uint64_t> st;
   vector<Graph::VertexPtr> ret;

   neighbors( vtx, depth, ret, st );
   ret.erase( ret.begin() );

   return ret;
}

void Graph::neighbors( const VertexPtr& vtx, int depth, vector<Graph::VertexPtr>& v, unordered_set<uint64_t>& st ) const
{
   uint64_t id = vtx.id();
   if ( !st.count( id ) )
   {
      st.insert( id );
      v.push_back( vtx );
   }

   if ( depth <= 0 )
      return;

   vector<Graph::VertexPtr> ret;

   for ( const Graph::VertexPtr& neighb : _Vertices[vtx._Index]._Neighbors )
      neighbors( neighb.premul( vtx._Mtx ), depth-1, v, st );
}

int Graph::colorOf( const TilePtr& tile ) const
{
   Perm perm = _Ico.colorPermOf( tile._Mtx );
   return perm[_Tiles[tile._Index]._Color];
}

vector<int> Graph::colorsAt( const VertexPtr& vtx ) const
{
   vector<int> ret;
   for ( const TilePtr& tile : tilesAt( vtx ) )
      ret.push_back( colorOf( tile ) );
   return ret;
}

uint32_t Graph::colorBits( const VertexPtr& vtx ) const
{
   return toBitFlag( colorsAt( vtx ) );
}

vector<Graph::TilePtr> Graph::tilesAt( const VertexPtr& vtx ) const
{
   if ( !vtx.isValid() )
      return {};

   vector<TilePtr> ret;   
   for ( const TilePtr& tile : _Vertices[vtx._Index]._Tiles )
      ret.push_back( tile.premul( vtx._Mtx ) );
   return ret;
}

Graph::TilePtr Graph::tileWithColor( const VertexPtr& vtx, int color ) const
{
   for ( const TilePtr& tile : tilesAt( vtx ) )
      if ( colorOf( tile ) == color )
         return tile;
   return TilePtr();
}

bool Graph::mustBeFar( const VertexPtr& a, const VertexPtr& b ) const
{
   for ( const TilePtr& tileA : tilesAt( a ) )
   {
      TilePtr tileB = tileWithColor( b, colorOf( tileA ) );
      if ( tileB.isValid() && !eq( tileA, tileB ) )
         return true;
   }
   return false;
}

bool Graph::mustBeClose( const VertexPtr& a, const VertexPtr& b ) const
{
   for ( const TilePtr& tileA : tilesAt( a ) )
   {
      TilePtr tileB = tileWithColor( b, colorOf( tileA ) );
      if ( tileB.isValid() && eq( tileA, tileB ) )
         return true;
   }
   return false;
}

bool Graph::eq( const TilePtr& a, const TilePtr& b ) const
{
   if ( a._Index != b._Index )
      return false;
   if ( _Tiles[a._Index]._SymmetryMap->match( a._Mtx, b._Mtx ) )
      return true;
   //return matrixId( a._Mtx ) == matrixId( b._Mtx );
   return false;
}

bool Graph::eq( const VertexPtr& a, const VertexPtr& b ) const
{
   if ( a._Index != b._Index )
      return false;
   if ( _Vertices[a._Index]._IsSymmetrical )
      return true;
   return matrixId( a._Mtx ) == matrixId( b._Mtx );
}

vector<Graph::VertexPtr> Graph::allVertices() const
{
   vector<VertexPtr> ret;
   for ( const IcoSymmetry::Config& config : _Ico.matrices() )
      for ( const Vertex& vtx : _Vertices )
         ret.push_back( VertexPtr( vtx._Index, config.m ) );
   return ret;
}

vector<Graph::VertexPtr> Graph::rawVertices() const
{
   vector<VertexPtr> ret;
   for ( const Vertex& vtx : _Vertices )
      ret.push_back( VertexPtr( vtx._Index, QMtx4x4() ) );
   return ret;
}

vector<Graph::TilePtr> Graph::rawTiles() const
{
   vector<TilePtr> ret;
   for ( int i = 0; i < (int) _Tiles.size(); i++ )
      ret.push_back( TilePtr( i, QMtx4x4() ) );
   return ret;
}

vector<Graph::TilePtr> Graph::allTiles() const
{
   vector<TilePtr> ret;
   for ( const IcoSymmetry::Config& config : _Ico.matrices() )
      for ( int i = 0; i < (int) _Tiles.size(); i++ )
         ret.push_back( TilePtr( i, QMtx4x4() ).premul( config.m ) );
   return ret;
}

vector<Graph::KeepCloseFar> Graph::calcKeepCloseFars() const
{
   vector<KeepCloseFar> ret;
   for ( const VertexPtr& vtx : rawVertices() )
   {      
      for ( const Graph::VertexPtr& neighb : neighbors( vtx, 5 ) )
      {
         KeepCloseFar kcf;
         kcf.a = vtx;
         kcf.b = neighb;
         kcf.keepClose = mustBeClose( vtx, neighb );
         kcf.keepFar = mustBeFar( vtx, neighb );
         if ( kcf.keepClose || kcf.keepFar )
            ret.push_back( kcf );
      }
   }
   return ret;
}

vector<Graph::LineVertexConstraint> Graph::calcLineVertexConstraints() const
{
   vector<LineVertexConstraint> ret;
   for ( const VertexPtr& a0 : rawVertices() )
   {      
      for ( const VertexPtr& a1 : neighbors( a0 ) ) if ( a0._Index <= a1._Index )
      {
         VertexPtr curveCenter = calcCurve( a0, a1 );
         for ( const TilePtr& tile : tilesAt( a0, a1 ) )
         {
            int color = colorOf( tile );
            for ( const Graph::VertexPtr& neighb : neighbors( a0, 6 ) )
            {
               TilePtr otherTile = tileWithColor( neighb, color );
               if ( !otherTile.isValid() )
                  continue;
               if ( eq( otherTile, tile ) )
                  continue;
               if ( eq( otherTile, tileWithColor( curveCenter, color ) ) )
                  continue; // otherTile is concave here
            
               LineVertexConstraint lvc;
               lvc.a0 = a0;
               lvc.a1 = a1;
               lvc.b = neighb;
               lvc.curveCenter = curveCenter;
               ret.push_back( lvc );
               
               if ( lvc.a0 == lvc.a1 )
                  qDebug() << "bad LineVertexConstraint";
            }
         }
      }
   }
   return ret;
}

vector<pair<Graph::VertexPtr, Graph::VertexPtr>> Graph::calcPerimeter() const
{
   vector<pair<VertexPtr, VertexPtr>> ret;
   
   //for ( const VertexPtr& vtx : rawVertices() )
   //{      
   //   for ( const Graph::VertexPtr& neighb : neighbors( vtx ) )
   //   {
   //      KeepCloseFar kcf;
   //      kcf.a = vtx;
   //      kcf.b = neighb;
   //      kcf.keepClose = !canBeFar( vtx, neighb );
   //      kcf.keepFar = !canBeClose( vtx, neighb );
   //      if ( kcf.keepClose || kcf.keepFar )
   //         ret.push_back( kcf );
   //   }
   //}
   for ( const Tile& tile : _Tiles )
   {
      for ( int i = 0; i < (int)tile._Vertices.size(); i++ )
      {
         VertexPtr a = tile._Vertices[i];
         VertexPtr b = tile._Vertices[(i+1)%tile._Vertices.size()];
         bool aOnPerim = false;
         bool bOnPerim = false;
         for ( const TilePtr& t : tilesAt( a ) ) if ( matrixId( t._Mtx ) != matrixId( QMtx4x4() ) ) aOnPerim = true;
         for ( const TilePtr& t : tilesAt( b ) ) if ( matrixId( t._Mtx ) != matrixId( QMtx4x4() ) ) bOnPerim = true;
         if ( aOnPerim && bOnPerim )
            ret.push_back( {a,b} );
      }
   }

   return ret;
}


vector<Graph::TilePtr> Graph::tilesAt( const VertexPtr& a, const VertexPtr& b ) const
{
   vector<TilePtr> ret;
   for ( const TilePtr& tile : tilesAt( a ) )
   {
      bool bIsAlsoOnTile = false;
      for ( const TilePtr& tileB : tilesAt( b ) )
         if ( eq( tile, tileB ) )
            bIsAlsoOnTile = true;
      if ( bIsAlsoOnTile )
         ret.push_back( tile );
   }
   return ret;
}

vector<Graph::VertexPtr> Graph::verticesForTile( const TilePtr& tile ) const
{
   vector<Graph::VertexPtr> ret;
   for ( const VertexPtr& vtx : _Tiles[tile._Index]._Vertices )
      ret.push_back( vtx.premul( tile._Mtx ) );
   return ret;
}

bool Graph::VertexPtr::operator==( const VertexPtr& rhs ) const
{
   return _Index == rhs._Index && matrixId( _Mtx ) == matrixId( rhs._Mtx );
}

Graph::VertexPtr Graph::calcCurve( const VertexPtr& a, const VertexPtr& b ) const
{
   vector<TilePtr> tiles = tilesAt( a, b );
   if ( tiles.size() != 2 )
      return VertexPtr();

   for ( int tileIdx = 0; tileIdx < 2; tileIdx++ )
   {
      int otherTileColor = colorOf( tiles[1-tileIdx] );
      for ( const VertexPtr& vtx : verticesForTile( tiles[tileIdx] ) )
      {
         if ( vtx == a ) continue;
         if ( vtx == b ) continue;
         //if ( mustBeFar( vtx, a ) && mustBeFar( vtx, b ) )
         //   return vtx; // wrong!
         vector<int> colors = colorsAt( vtx );
         if ( find( colors.begin(), colors.end(), otherTileColor ) != colors.end() )
            return vtx;
      }
   }
   return VertexPtr();
}







bool Dual::Vertex::hasNeighbor( const VertexPtr& a ) const
{
   return findNeighborIdx( a ) >= 0;
}

void Dual::Vertex::eraseNeighbor( const VertexPtr& a )
{
   int idx = findNeighborIdx( a );
   if ( idx < 0 )
      return;
   _Neighbors.erase( _Neighbors.begin() + idx );
}

int Dual::Vertex::findNeighborIdx( const VertexPtr& a ) const
{
   for ( int i = 0; i < (int)_Neighbors.size(); i++ )
      if ( _Neighbors[i] == a )
         return i;
   return -1;
}

bool Dual::VertexPtr::operator==( const VertexPtr& rhs ) const
{
   return _Index == rhs._Index && matrixId( _Mtx ) == matrixId( rhs._Mtx );
}


bool Dual::isDuplicate( const VertexPtr& a ) const
{
   return !_Vertices[a._Index]._SymmetryMap->isReal( a._Mtx );
}

Dual::VertexPtr Dual::toReal( const VertexPtr& a ) const
{
   return VertexPtr( a._Index, _Vertices[a._Index]._SymmetryMap->toReal( a._Mtx ) );
}

vector<Dual::VertexPtr> Dual::allVertices() const
{
   vector<Dual::VertexPtr> ret;
   for ( const Vertex& vertex : _Vertices )
   for ( const IcoSymmetry::Config& config : _Ico.matrices() )
   {
      VertexPtr a( vertex._Index, config.m );
      if ( !isDuplicate( a ) )
         ret.push_back( a );
   }
   return ret;
}

XYZ Dual::posOf( const VertexPtr& vtx ) const
{
   return vtx._Mtx * _Vertices[vtx._Index]._Pos;
}


void Dual::setPos( const VertexPtr& vtx, const XYZ& pos )
{
   _Vertices[vtx._Index]._Pos = vtx._Mtx.inverted() * pos;
   //return vtx._Mtx * _Vertices[vtx._Index]._Pos;
}

void Dual::toggleEdge( const VertexPtr& a, const VertexPtr& b, bool onlyAdd )
{
   if ( !a.isValid() || !b.isValid() )
      return;
   
   if ( a == b )
   {
      qDebug( "toggleEdge dup" );
      return;
   }

   VertexPtr bb = toReal( premul( b, a._Mtx.inverted() ) );
   VertexPtr aa = toReal( premul( a, b._Mtx.inverted() ) );

   if ( _Vertices[a._Index].hasNeighbor( bb ) && _Vertices[b._Index].hasNeighbor( aa ) )
   {
      if ( !onlyAdd )
      {
         _Vertices[a._Index].eraseNeighbor( bb );
         _Vertices[b._Index].eraseNeighbor( aa );
      }
   }
   else
   {   
      _Vertices[a._Index]._Neighbors.push_back( bb );
      if ( idOf( aa ) != idOf( bb ) ) // don't double-add symmetric edge
         _Vertices[b._Index]._Neighbors.push_back( aa ); 
   }
}

vector<Dual::VertexPtr> Dual::baseVertices() const
{
   vector<VertexPtr> ret;
   for ( const Vertex& a : _Vertices )
      ret.push_back( a.toVertexPtr() );
   return ret;
}

vector<Dual::VertexPtr> Dual::neighborsOf( const VertexPtr& a ) const
{
   if ( !a.isValid() )
      return {};

   vector<VertexPtr> ret;
   for ( const QMtx4x4& mtx : _Vertices[a._Index]._SymmetryMap->symmetricMatrices( a._Mtx ) )
      for ( const VertexPtr& b : _Vertices[a._Index]._Neighbors )
         ret.push_back( premul( b, mtx ) );
   return ret;
}

vector<Dual::VertexPtr> Dual::sortedNeighborsOf( const VertexPtr& a ) const
{
   vector<Dual::VertexPtr> v = neighborsOf( a );
   
   QMtx4x4 m = matrixRotateToZAxis( posOf( a ) );
   auto angleOf = [&]( const XYZ& p ) { XYZ q = m*p; return ::atan2( q.y, q.x ); };
   
   sort( v.begin(), v.end(), [&]( const Dual::VertexPtr& a, const Dual::VertexPtr& b ) { return angleOf( posOf( a ) ) < angleOf( posOf( b ) ); } );
   return v;
}

int Dual::colorOf( const VertexPtr& vtx ) const
{
   Perm perm = _Ico.colorPermOf( vtx._Mtx );
   return perm[_Vertices[vtx._Index]._Color];
}

void Dual::setColorOf( const VertexPtr& vtx, int color )
{
   Perm perm = _Ico.colorPermOf( vtx._Mtx );
   _Vertices[vtx._Index]._Color = perm.inverted()[color];
}

int Dual::idOf( const VertexPtr& a ) const
{
   VertexPtr aa = toReal( a );
   return MatrixIndexMap::indexOf( a._Mtx ) * _Vertices.size() + aa._Index;
}

Dual::VertexPtr Dual::fromId( int id ) const
{
   int sz = (int) _Vertices.size();   
   return VertexPtr( id % sz, MatrixIndexMap::at( id / sz ) );
}

Dual::VertexPtr Dual::next( const VertexPtr& a, const VertexPtr& b ) const
{
   vector<VertexPtr> v = sortedNeighborsOf( a );
   for ( int i = 0; i < (int)v.size(); i++ )
      if ( idOf( v[i] ) == idOf( b ) )
         return v[(i+1)%(int)v.size()];
   return VertexPtr();
}

vector<Dual::VertexPtr> Dual::polygon( const VertexPtr& a, const VertexPtr& b ) const
{
   vector<VertexPtr> ret = { b, a };

   while ( true )
   {
      VertexPtr c = next( ret.back(), ret[ret.size()-2] );
      if ( idOf( c ) == idOf( ret[0] ) )
         return ret;
      ret.push_back( c );
   }
}

Dual::VertexPtr Dual::premul( const VertexPtr& vtx, const QMtx4x4& mtx ) const 
{ 
   return toReal( VertexPtr( vtx._Index, mtx * vtx._Mtx ) ); 
}