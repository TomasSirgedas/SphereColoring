#include "Drawing.h"
#include <QPainter>
#include <QMatrix4x4>
#include <QBrush>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <QMouseEvent>

#include "Model.h"
#include "Simulation.h"

using namespace std;

const vector<uint32_t> COLORS = { 0xFFED1D26, 0xFF3F47CB, 0xFFFFF200, 0xFFA349A4, 0xFF22B14C, 0xFFFF7F27, 0xFF00EEEE, 0xFFFFFFFF, 0xFF000000, 0xFF808080 };

//int debugPermute[60][8] = {0};
//vector<vector<int>> debugPermute = {{0, 0, 0, 0, 0, 0, 0, 0}, {6, 7, 3, 6, 6, 6, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {6, 4, 6, 6, 0, 5, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {3, 4, 0, 5, 2, 2, 5, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {5, 5, 5, 6, 0, 5, 2, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {9, 0, 3, 3, 3, 6, 4, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {6, 4, 6, 0, 5, 6, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {6, 3, 13, 6, 6, 0, 1, 0}, {8, 3, 4, 4, 6, 7, 10, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 5, 0, 0, 1, 2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {4, 6, 6, 0, 5, 8, 6, 0}, {6, 4, 0, 5, 6, 6, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}};
//vector<vector<int>> debugPermute( 60, vector<int>( 8 ));

//Perm finalColorPerm({6, 2, 0, 3, 4, 1, 5});
//QColor tileColor( int idx ) { vector<int> perm = {5, 1, 0, 4, 3, 2, 6}; return QColor::fromRgba( COLORS[perm[idx]] ); }
QColor tileColor( int idx ) { return QColor::fromRgba( COLORS[idx] ); }
QColor withAlpha( const QColor& color, int alpha ) { return QColor( color.red(), color.green(), color.blue(), alpha ); }

const double PI = acos(0.) * 2.;


XYZ lineSphereIntersection( const XYZ& p0, const XYZ& p1, double radius )
{
   XYZ v = p1-p0;
   double a = v*v;
   double b = p0*v;
   double c = p0*p0 - radius*radius;
   double t = ( -b - sqrt( b*b - a*c ) ) / a;
   return p0 + v * t;
}

//vector<XYZ> calcCurve( const XYZ& a, const XYZ& b, const XYZ& center, double maxDistance )
//{
//   int numSegments = (int) ceil( a.dist( b ) / maxDistance );
//   if ( numSegments <= 1 )
//      return { b };
//   
//   XYZ mid = ((a-center)+(b-center)).normalized() * (( (a-center).len() + (b-center).len() )/2) + center;   
//   vector<XYZ> part0 = calcCurve( a, mid, center, maxDistance );
//   vector<XYZ> part1 = calcCurve( mid, b, center, maxDistance );
//   vector<XYZ> ret = part0;
//   ret.insert( ret.end(), part1.begin(), part1.end() );
//   return ret;
//}

// dir == 1, dir == -1 --> curve direction
// dir == 0 --> pick shorter curve direction
vector<XYZ> calcCurve2( const XYZ& p0, const XYZ& p1, const XYZ& center, double maxDistance, int dir )
{
   if ( p0.dist2(p1) < 1e-14 )
      return {};
   XYZ a = center.len2() < 1e-14 ? (p0^p1).normalized() : center.normalized();
   if ( (a^p0).len2() < 1e-14 )
      return {}; // p is on the axis
      
   XYZ v = (a^p0).normalized();
   XYZ u = v^a;
   double angle = atan2( p1*v, p1*u );
   if ( dir == 0 ) dir = angle > 0 ? 1 : -1;
   angle *= dir;
   if ( angle < 0 )
      angle += PI*2;

   double radius0 = p0 * u;
   double radius1 = sqrt((p1*u)*(p1*u) + (p1*v)*(p1*v));
   double dist = radius0 * angle;
   int numSegments = (int) ceil( dist / maxDistance );

   vector<XYZ> ret;
   for ( int i = 1; i <= numSegments; i++ )
   {
      double t = (double)i / numSegments;
      double radius = radius0 * (1-t) + radius1 * t;
      double zDist = (p0*a)*(1-t) + (p1*a)*t;
      ret.push_back( a*zDist + u*radius*cos(dir*angle*t) + v*radius*sin(dir*angle*t) );
   }
   return ret;
}

vector<XYZ> calcPolyCurveOnSphere( const vector<XYZ>& v, double maxDistance, int dir )
{
   vector<XYZ> ret;
   for ( int i = 0; i < (int) v.size(); i++ )
   {
      vector<XYZ> segment = calcCurve2( v[i], v[(i+1)%v.size()], XYZ(0,0,0), maxDistance, dir );
      ret.insert( ret.end(), segment.begin(), segment.end() );
   }
   return ret;
}

// returns point at a distance 1 from a and b, on the sphere
XYZ calcOutlinePt( const XYZ& a, const XYZ& b, double radius )
{
   double u = (radius*radius-.5) / (radius*radius + a*b);
   double v = sqrt( (radius*radius - (a+b)*(a+b)*u*u) / ((a^b)*(a^b)) );
   return (a+b)*u + (a^b)*v;
}

XYZ rotateOnAxis( const XYZ& p, const XYZ& axis, double angle )
{
   XYZ a = axis.normalized();
   if ( (a^p).len2() < 1e-14 )
      return p; // p is on the axis
   XYZ v = (a^p).normalized();
   XYZ u = a^v;

   return a*(p*a) + u*(p*u)*cos(angle) + v*(p*u)*sin(angle);
}

vector<XYZ> calcTileOutline( const Graph& graph, const Graph::TilePtr& tile, double maxSpacing )
{
   vector<XYZ> ret;

   vector<Graph::VertexPtr> v = graph.verticesForTile( tile );
   for ( int i = 0; i < (int)v.size(); i++ )
   {
      if ( maxSpacing >= 1 ) { ret.push_back( graph.posOf( v[i] ) ); continue; }

      const Graph::VertexPtr& a = v[i];
      const Graph::VertexPtr& b = v[(i+1)%v.size()];
      Graph::VertexPtr c = graph.calcCurve( a, b ); // c = center of curve
      vector<XYZ> curve;
      if ( c.isValid() )
         curve = calcCurve2( graph.posOf( a ), graph.posOf( b ), graph.posOf( c ), maxSpacing, 0 );
      else
         curve = calcCurve2( graph.posOf( a ), graph.posOf( b ), XYZ(), maxSpacing, 1 );
      ret.insert( ret.end(), curve.begin(), curve.end() );
   }

   return ret;
}

// dir == 1, dir == -1 --> curve direction
// dir == 0 --> pick shorter curve direction
vector<XYZ> calcCurve( const XYZ& p0, const XYZ& p1, const XYZ& center, double maxDistance, int dir )
{
   if ( p0.dist2(p1) < 1e-14 )
      return {};
   XYZ a = center.len2() < 1e-14 ? (p0^p1).normalized() : center.normalized();
   if ( (a^p0).len2() < 1e-14 )
      return {}; // p is on the axis

   XYZ v = (a^p0).normalized();
   XYZ u = v^a;
   double angle = atan2( p1*v, p1*u );
   if ( dir == 0 ) dir = angle > 0 ? 1 : -1;
   angle *= dir;
   if ( angle < 0 )
      angle += PI*2;

   double radius0 = p0 * u;
   double radius1 = sqrt((p1*u)*(p1*u) + (p1*v)*(p1*v));
   double dist = radius0 * angle;
   int numSegments = (int) ceil( dist / maxDistance );

   vector<XYZ> ret;
   for ( int i = 1; i <= numSegments; i++ )
   {
      double t = (double)i / numSegments;
      double radius = radius0 * (1-t) + radius1 * t;
      double zDist = (p0*a)*(1-t) + (p1*a)*t;
      ret.push_back( a*zDist + u*radius*cos(dir*angle*t) + v*radius*sin(dir*angle*t) );
   }
   return ret;
}

vector<XYZ> expandOutlineOnSphere( const vector<XYZ>& v, double maxSpacing )
{   
   if ( v.empty() ) 
      return {};
   vector<XYZ> ret;

   double R = v[0].len();
   
   for ( int i0 = 0; i0 < (int) v.size(); i0++ )
   {
      int i1 = (i0+1) % v.size();
      int i2 = (i0+2) % v.size();

      XYZ p0 = calcOutlinePt( v[i1], v[i0], R );
      XYZ p1 = calcOutlinePt( v[i2], v[i1], R );

      vector<XYZ> curve = calcCurve2( p0, p1, v[i1], maxSpacing, 0 );
      ret.insert( ret.end(), curve.begin(), curve.end() );
   }

   return ret;
}

QPolygonF toQPolygonF( vector<XYZ> path, const QMtx4x4& modelToBitmap )
{   
   //QMtx4x4 rotatedModelToBitmap = calcRotatedModeltoBitmap( r() );
   function<QPointF(const XYZ&)> rotatedToBitmap = [&]( const XYZ& p ) { XYZ bitmapPos = modelToBitmap * p; return QPointF( bitmapPos.x, bitmapPos.y ); };

   QPolygonF poly;
   XYZ edgePt;

   // rotate path so that tilePath.end() -> tilePath.begin() is visible
   {
      int startIdx = -1;
      for ( int i = 0; i < (int) path.size(); i++ )
         if ( path[i].z < 0 && path[(i+1)%path.size()].z < 0  )
         { startIdx = (i+1)%path.size(); break; }
      if ( startIdx < 0 )
         return QPolygonF();
      rotate( path.begin(), path.begin() + startIdx, path.end() );
   }

   for ( const XYZ& p : path )
   {
      if ( p.z < 0 )
      {
         if ( edgePt != XYZ() )
         {
            for ( const XYZ& q : calcCurve( edgePt, XYZ( p.x, p.y, 0 ).normalized() * p.len(), XYZ(), .1, 1 ) )
               poly.append( rotatedToBitmap( q ) );                     
            edgePt = XYZ();
         }

         poly.append( rotatedToBitmap( p ) );
      }
      else
      {
         if ( edgePt == XYZ() )
            edgePt = XYZ( p.x, p.y, 0 ).normalized() * p.len();
      }
   }

   return poly;
}


Drawing::Drawing( QWidget *parent )
   : QWidget( parent )
{
   ui.setupUi( this );
}

Drawing::~Drawing()
{
}


void Drawing::resizeEvent( QResizeEvent *event )
{
   updateLabel();
}

void Drawing::updateLabel()
{         
   if ( _Simulation )
      ui.label->setPixmap( QPixmap::fromImage( makeImage( *_Simulation->_Graph ) ) );
}

//void Drawing::mousePressEvent( QMouseEvent* event )
//{
//   _ClickPos = event->pos();
//   updateLabel();
//}

QMtx4x4 Drawing::modelToBitmap() const
{
   return modelToBitmapNoRot() * modelRotation();
}

QMtx4x4 Drawing::modelToBitmapNoRot() const
{
   QSize size = this->size();
   QMtx4x4 ret;
   ret.translate( size.width() / 2, size.height() / 2, 0 );
   ret.scale( size.height() / _Radius * .49 );
   ret.scale( 1, -1 );
   ret.scale( _Zoom );
   return ret;
}

QMtx4x4 Drawing::modelRotation() const
{
   QMtx4x4 ret;
   return _ModelRotation;
   //ret.rotateX( -_XRotation/60. );
   //ret.rotateY( _YRotation/60. );
   //return ret;
}


bool Drawing::isOnNearSide( const XYZ& p ) const
{
   QMtx4x4 modelToBitmap = Drawing::modelToBitmap();
   XYZ a = modelToBitmap * p;
   return a.z < 0;
}

bool Drawing::bitmapToModel( const QPoint& p, XYZ& modelPos ) const
{
   QMtx4x4 modelToBitmap = Drawing::modelToBitmap();
   XYZ p0 = modelToBitmap.inverted() * XYZ( p.x(), p.y(), 0 );
   XYZ p1 = modelToBitmap.inverted() * XYZ( p.x(), p.y(), 1 );
   modelPos = lineSphereIntersection( p0, p1, _Radius );
   return true;
}

QImage Drawing::makeImage( Graph& graph )
{
   if ( !_ShowDual && &graph == nullptr )
      return QImage();

   double radius = _Radius;
   QSize size = this->size();
   QImage image( size, QImage::Format_RGB888 );
   image.fill( QColor( 0, 0, 0 ) );

   shared_ptr<Dual> dual = _Simulation->_Dual;

   QPainter painter( &image );

   painter.setRenderHint( QPainter::Antialiasing, true );   

   QMtx4x4 modelToBitmap = Drawing::modelToBitmap();

   //IcoSymmetry ico;
   //XYZ ico0 = ico[0].normalized() * radius;
   //XYZ ico01 = (( ico[0] + ico[1] ) / 2).normalized() * radius;
   //XYZ ico02 = (( ico[0] + ico[2] ) / 2).normalized() * radius;
   //XYZ ico012 = (( ico[0] + ico[1] + ico[2] ) / 3).normalized() * radius;
   //XYZ ico015 = (( ico[0] + ico[1] + ico[5] ) / 3).normalized() * radius;
   //XYZ ico1 = ico[1].normalized() * radius;



   //TileDots tileDots = generateSphereColoringDots( _Custom[0] );
   //Graph graph( tileDots );
   //int EXTENSIONS = _Custom[0]-1;
   //Graph graph( EXTENSIONS );


   QPoint clickedPos = _ClickPos;
   Graph::VertexPtr clickedVtx;
   if ( _ClickPos.x() > 0 )
   {
      Graph::VertexPtr bestVtx;
      double bestDist = 9999;
      for ( const ISymmetry::Config& config : GlobalSymmetry::matrices() )
      for ( const Graph::VertexPtr& a_ : graph.rawVertices() )
      {
         Graph::VertexPtr a( a_._Index, config.m );
         if ( graph.posOf( a ).z >= 0 )
            continue;
         QPointF bitmapPos = ( modelToBitmap * graph.posOf( a ) ).toPointF();
         double dist = QLineF( bitmapPos, _ClickPos ).length();
         if ( dist >= bestDist )
            continue;
         bestDist = dist;
         bestVtx = a;
      }
      clickedVtx = bestVtx;
      _ClickPos = QPoint(0,0);
      //qDebug() << graph._Vertices[clickedVtx._Index]._HexPos;
   }

   if ( clickedVtx.isValid() )
   {
      XYZ p;
      bitmapToModel( clickedPos, p );

      graph._Vertices[clickedVtx._Index]._Pos = clickedVtx._Mtx.inverted() * p;
   }

   if ( clickedVtx.isValid() )
   {
      _SelectedVertices.push_back( clickedVtx );
      while ( _SelectedVertices.size() > 2 )
         _SelectedVertices.erase( _SelectedVertices.begin() );
   }


   for ( int stage : { 1, 3, 4, 6, 7, 8, 9, 10, 11 } )
      for ( const ISymmetry::Config& config : GlobalSymmetry::matrices() )
      {
         const QMtx4x4& m = config.m;

         auto toBitmap = [&]( const XYZ& pos ) { return ( modelToBitmap * m * pos ).toPointF(); };
         auto toBitmapNoRotate = [&]( const XYZ& pos ) { return ( modelToBitmap * pos ).toPointF(); };
                  
         if ( stage == 1 && !_ShowDual )
         {
            Perm quadPerm = GlobalSymmetry::colorPermOf( m );

            int alpha = 255;
            if ( config.isHomeState() )
               alpha = 255;

            int idx = 0;
            for ( const Graph::TilePtr& tile_ : graph.rawTiles() )
            {
               Graph::TilePtr tile = tile_.premul( m );

               bool tileVisible = false;
               for ( const Graph::VertexPtr& a : graph.verticesForTile( tile ) )
                  if ( isOnNearSide( graph.posOf( a ) ) )
                     tileVisible = true;
               if ( !tileVisible )
                  continue;

               painter.setPen( Qt::NoPen );
               painter.setBrush( withAlpha( tileColor( graph.colorOf(tile) ), alpha ) );

               //{
               //   QPolygonF poly;                  
               //   for ( const Graph::VertexPtr& vtx : tile._Vertices )
               //      poly.append( toBitmap( graph.posOf( vtx ) ) );
               //   painter.drawPolygon( poly );
               //}
               {
                  QPolygonF poly;                  
                  for ( const XYZ& p : calcTileOutline( graph, tile, _DrawCurves ? .1 : 1 ) )
                     poly.append( toBitmapNoRotate( p ) );
                  painter.drawPolygon( poly );
               }
               idx++;
            }
         }
         if ( stage == 2 && !_ShowDual )
         {
            //painter.setPen( QPen( QColor(255,255,255,128), 3 ) );
            painter.setPen( QPen( QColor(0,0,0,96), 3 ) );
            painter.setBrush( Qt::NoBrush );
            for ( const auto& pr : _Simulation->_Graph->calcPerimeter() )
            {
               XYZ a = _Simulation->_Graph->posOf( pr.first );
               XYZ b = _Simulation->_Graph->posOf( pr.second );
               painter.drawLine( toBitmap( a ), toBitmap( b ) );
            }
         }
         if ( stage == 3 && config.isHomeState() && _ShowViolations && !_ShowDual )
         {
            painter.setBrush( QColor( 0, 0, 0, 64 ) );
            painter.setPen( Qt::NoPen );
            painter.drawRect( image.rect() );

            const double TOLERANCE = 1e-5;

            painter.setBrush( Qt::NoBrush );
            for ( const ISymmetry::Config& config : GlobalSymmetry::matrices() )
            for ( const auto& pr : _Simulation->_KeepCloseFars )
            {
               XYZ a = config.m * _Simulation->_Graph->posOf( pr.a );
               XYZ b = config.m * _Simulation->_Graph->posOf( pr.b );
               if ( !isOnNearSide( a ) || !isOnNearSide( b ) )
                  continue;
               double dist = a.dist( b );
               if ( pr.keepFar && dist < 1. - TOLERANCE )
               {
                  painter.setPen( QPen( QColor(0,255,0,96), 3 ) );
                  painter.drawLine( toBitmap( a ), toBitmap( b ) );
               }
               if ( pr.keepClose && dist > 1. + TOLERANCE )
               {
                  painter.setPen( QPen( QColor(255,0,0,96), 3 ) );
                  painter.drawLine( toBitmap( a ), toBitmap( b ) );
               }
            }
         }    
         if ( stage == 4 && config.isHomeState() && _DrawRigidEDs && !_ShowDual )
         {
            painter.setPen( QPen( QColor(0,0,0,128), 1 ) );
            painter.setBrush( Qt::NoBrush );
            for ( const ISymmetry::Config& config : GlobalSymmetry::matrices() )
               for ( const auto& pr : _Simulation->_KeepCloseFars )
               {  
                  XYZ a = graph.posOf( pr.a.premul( config.m ) );
                  XYZ b = graph.posOf( pr.b.premul( config.m ) );
                  if ( pr.keepClose && pr.keepFar && isOnNearSide( a ) && isOnNearSide( b )  )
                     painter.drawLine( toBitmap( a ), toBitmap( b ) );
               }
         }     
         if ( stage == 5 && config.isHomeState() && _LabelVertices && !_ShowDual )
         {            
            painter.setPen( QColor( 0, 0, 0, 255 ) );
            for ( const Graph::VertexPtr& a : graph.rawVertices() )
            {
               painter.drawText( toBitmap( graph.posOf( a ) ), QString::number( a._Index ) );               
            }
         }
         if ( stage == 6 && config.isHomeState() && _DrawZoneOfExclusions && !_ShowDual )
         {
            painter.setPen( QPen( QColor(255,255,255,192), 0 ) );
            painter.setBrush( Qt::NoBrush );

            for ( const Graph::TilePtr& tile : graph.allTiles() ) if ( graph.colorOf( tile ) == lround(_Custom[0]) )
            //const Graph::Tile& tile = graph._Tiles[lround(_Custom[0])];
            {
               if ( !graph._Tiles[tile._Index]._SymmetryMap->isReal( tile._Mtx ) ) // only use one copy
                  continue;

               vector<XYZ> outline = expandOutlineOnSphere( calcTileOutline( graph, tile, .02 ), .02 );
               bool allOnNearSide = true;
               for ( const XYZ& p : outline )
                  if ( !isOnNearSide( p ) )
                     allOnNearSide = false;
               if ( !allOnNearSide )
                  continue;

               QPolygonF poly;
               for ( const XYZ& p : outline )
                  poly.append( toBitmap( p ) );               
               painter.drawPolygon( poly );

               //for ( const XYZ& p : expandOutlineOnSphere( calcTileOutline( graph, tile, .02 ), .02 ) )
               //   painter.drawEllipse( toBitmap( p ), 2, 2 );
               //for ( const XYZ& p : calcTileOutline( graph, tile, .02 ) )
               //   painter.drawEllipse( toBitmap( p ), 2, 2 );
            }
         }
         // draw outline of 1/60th sector
         if ( stage == 7 && _DrawSectors )
         {
            painter.setPen( QPen( QColor(255,255,255,192), 1 ) );
            painter.setBrush( Qt::NoBrush );
            vector<XYZ> polyCurve = calcPolyCurveOnSphere( GlobalSymmetry::sectorOutline( radius ), .1, 1 );
            painter.drawPolygon( toQPolygonF( modelRotation() * m * polyCurve, modelToBitmapNoRot() ) );
         }
         // edges
         if ( stage == 9 && config.isHomeState() && _ShowDual )
         {
            painter.setBrush( Qt::NoBrush );

            //for ( const GlobalSymmetry::Config& config : GlobalSymmetry::matrices() )
            {
               for ( const Dual::VertexPtr& a : dual->allVertices() )
               {
                  painter.setPen( QPen( QColor(255,255,255, 16), 1 ) );
                  for ( const Dual::VertexPtr& b : dual->neighborsOf( a ) )
                  {
                     XYZ posA = config.m * dual->posOf( a );
                     XYZ posB = config.m * dual->posOf( b );
                     if ( isOnNearSide( posA ) && isOnNearSide( posB ) )
                        painter.drawLine( toBitmap( posA ), toBitmap( posB ) );
                  }
               }
            }
         }
         // vertices
         if ( stage == 8 && config.isHomeState() && _ShowDual )
         {
            painter.setBrush( Qt::NoBrush );

            for ( const Dual::VertexPtr& a : dual->allVertices() )
            {
               int color = dual->colorOf( a );
               //color = (color + debugPermute[MatrixIndexMap::indexOf(a._Mtx)][color]) % 7;
               painter.setBrush( tileColor( color ) );
               XYZ pos = dual->posOf( a );
               if ( isOnNearSide( pos ) )
               {
                  painter.setPen( Qt::black );
                  double r = dual->_Vertices[a._Index]._SymmetryMap->_SymmetricMatrices[0].size() == 5 ? 10 : 4;
                  painter.drawEllipse( toBitmap( pos ), r, r );
                  painter.setPen( QColor( 255, 255, 255, 64 ) );
                  if ( _LabelVertices )
                     painter.drawText( toBitmap( pos ) + QPointF( 0, -4 ), QString::number( dual->idOf( a ) ) );


                  //if ( dual->_Vertices[a._Index]._SymmetryMap->_SymmetricMatrices[0].size() == 5 )
                  //   painter.drawText( toBitmap( pos ) + QPointF( 0, -12 ), QString::number( ico.id( dual->posOf( a ).normalized() * 1.90211303259 ) ) );
               }
            }
         }
         // tile vertices
         if ( stage == 10 && config.isHomeState() && _LabelVertices && !_ShowDual )
         {
            painter.setBrush( Qt::NoBrush );
            painter.setPen( QColor( 0,0,0 ) );

            for ( const Graph::VertexPtr& a : graph.allVertices() ) if ( graph.posOf( a ).z < 0 )
            {
               //painter.drawEllipse( toBitmap( graph.posOf( a ) ), 2, 2 );
               painter.drawText( toBitmap( graph.posOf( a ) ) + QPointF( 0, 0 ), QString::number( graph.idOf( a ) ) );
            }
         }
         //// draw sector text
         //if ( stage == 11 && config.isHomeState() )
         //{
         //   painter.setBrush( Qt::NoBrush );
         //   painter.setPen( QColor( 255,255,255 ) );

         //   for ( const GlobalSymmetry::Config& config : GlobalSymmetry::matrices() )
         //   {
         //      XYZ p = config.m * ( ico[0]*3 + ico[1] + ico[2] ).normalized() * _Radius;
         //      if ( !isOnNearSide( p ) )
         //         continue;
         //      painter.drawText( toBitmap( p ), QString::number(MatrixIndexMap::indexOf(config.m)) );
         //   }

         //   //for ( const GlobalSymmetry::Config& config : GlobalSymmetry::matrices() )
         //   //   for ( const auto& pr : _Simulation->_KeepCloseFars )
         //   //   {  
         //   //      XYZ a = graph.posOf( pr.a.premul( config.m ) );
         //   //      XYZ b = graph.posOf( pr.b.premul( config.m ) );
         //   //      if ( pr.keepClose && pr.keepFar && isOnNearSide( a ) && isOnNearSide( b )  )
         //   //         painter.drawLine( toBitmap( a ), toBitmap( b ) );
         //   //   }
         //}
      }

   painter.setFont( QFont( "Arial", 24 ) );
   painter.setPen( QColor( 200, 200, 200 ) );
   painter.drawText( 60, image.height() - 15, "r = " + QString::number( _Simulation->_Radius ) );


   return image;
}

QPoint Drawing::mousePos() const
{
   return mapFromGlobal( QCursor::pos() );
}

Dual::VertexPtr Drawing::dualVertexNearest( const QPoint& mousePos, double maxPixelDist )
{
   double bestDist = maxPixelDist;
   Dual::VertexPtr bestVtx;
   
   QMtx4x4 modelToBitmap = Drawing::modelToBitmap();
   
   for ( const Dual::VertexPtr& vtx : _Simulation->_Dual->allVertices() )
   {
      XYZ p = modelToBitmap * _Simulation->_Dual->posOf( vtx );
      if ( p.z > 0 )
         continue;
      QPointF vertexPos = p.toPointF();
      double dist = QLineF( vertexPos, mousePos ).length();
      if ( dist < bestDist )
      {
         bestDist = dist;
         bestVtx = vtx;
      }
   }

   return bestVtx;
}

void Drawing::debugClick( QMouseEvent* event )
{
   //if ( !event->buttons().testFlag( Qt::LeftButton ) )
   //   return;

   //Dual::VertexPtr vtx = dualVertexNearest( event->pos(), 8 );
   //if ( !vtx.isValid() )
   //   return;
   ////qDebug() << _Simulation->_Dual->idOf( vtx ) << " " << MatrixIndexMap::indexOf( vtx._Mtx );
   //debugPermute[MatrixIndexMap::indexOf( vtx._Mtx )][_Simulation->_Dual->colorOf( vtx )]++;
   //qDebug() << debugPermute;
   //updateLabel();
}