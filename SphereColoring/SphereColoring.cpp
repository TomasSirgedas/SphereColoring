#include "SphereColoring.h"
#include "Drawing.h"
#include "Model.h"
#include "PlatformSpecific.h"
#include <QDebug>
#include <QShortcut>
#include <QMouseEvent>
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <QFile>
#include <QFileDialog>

#include <vector>
#include <sstream>
#include <iomanip>


using namespace std;


shared_ptr<Graph> makeGraph( shared_ptr<const Dual> dual, double radius )
{
   shared_ptr<Graph> graph( new Graph );
   
   std::map<set<int>, int> polygonToTileIndex;

   for ( int k = 0; k < (int) dual->_Vertices.size(); k++ )
   {
      Dual::VertexPtr a = dual->fromId( k );

      Graph::Tile tile;
      tile._Index = (int) graph->_Tiles.size();
      tile._Color = dual->colorOf( a );
      tile._SymmetryMap = dual->_Vertices[a._Index]._SymmetryMap;

      for ( const Dual::VertexPtr& b : dual->sortedNeighborsOf( a ) )
      {
         vector<Dual::VertexPtr> poly = dual->polygon( a, b );

         Graph::VertexPtr tileVertex;
         set<int> polyAsSet;
         for ( const Dual::VertexPtr& c : poly )
            polyAsSet.insert( dual->idOf( c ) );

         {
            for ( const IcoSymmetry::Config& config : GlobalSymmetry::matrices() )
            {
               set<int> polyAsSet;
               for ( const Dual::VertexPtr& c : poly )
                  polyAsSet.insert( dual->idOf( dual->premul( c, config.m ) ) );
               if ( polygonToTileIndex.count( polyAsSet ) )
               {
                  tileVertex = Graph::VertexPtr( polygonToTileIndex.at(polyAsSet), config.m.inverted() );
                  break;
               }
            }
         }         

         if ( !tileVertex.isValid() ) // create it if needed
         {
            XYZ sum;
            for ( const Dual::VertexPtr& c : poly )
               sum += dual->posOf( c );

            Graph::Vertex v( (int) graph->_Vertices.size() );
            v._IsSymmetrical = MatrixSymmetryMap::symmetryFor( sum )->hasSymmetry();
            v._Neighbors;
            v._Pos = sum.normalized() * radius;
            v._Tiles;
            graph->_Vertices.push_back( v );
            tileVertex = Graph::VertexPtr( v._Index, QMtx4x4() );
            polygonToTileIndex[polyAsSet] = v._Index;
         }

         tile._Vertices.push_back( tileVertex );         
      }
      graph->_Tiles.push_back( tile );
      for ( const Graph::VertexPtr& vtx : tile._Vertices )
         if ( tile._SymmetryMap->isReal( vtx._Mtx ) ) // only add one copy
            graph->_Vertices[vtx._Index]._Tiles.push_back( Graph::TilePtr( tile._Index, vtx._Mtx.inverted() ) );
      for ( int i = 0; i < (int) tile._Vertices.size(); i++ )
         graph->addNeighbor( tile._Vertices[i], tile._Vertices[(i+1)%tile._Vertices.size()] );
   }


//   for ( const Graph::VertexPtr& a : graph->allVertices() )
//      for ( const Graph::VertexPtr& b : graph->neighbors( a ) )
//         qDebug() << graph->idOf( a ) << "-" << graph->idOf( b );

   return graph;
}


void saveDual( const QString& filename, const Dual& dual )
{
   if ( filename.isEmpty() )
      return;
   
   QJsonArray vertices;
   for ( const Dual::Vertex& a : dual._Vertices )
      vertices.push_back( QJsonObject { { "color", a._Color }, { "x", a._Pos.x }, { "y", a._Pos.y }, { "z", a._Pos.z } } );

   QJsonArray edges;
   for ( const Dual::Vertex& a : dual._Vertices )
      for ( const Dual::VertexPtr& b : a._Neighbors ) if ( a._Index <= b._Index )
         edges.push_back( QJsonArray { a._Index, b._Index, MatrixIndexMap::indexOf( b._Mtx ) } );
   QJsonObject graph = { { "vertices", vertices }, { "edges", edges }, { "symmetry", QString::fromStdString( GlobalSymmetry::symmetry()->name() ) } };         
   {
      QFile f( filename );
      f.open(QFile::WriteOnly);
      f.write(QJsonDocument( graph ).toJson());            
   }   
}

shared_ptr<Dual> loadDual( const QString& filename )
{
   if ( filename.isEmpty() )
      return nullptr;

   shared_ptr<Dual> dual( new Dual );
   QFile f( filename );
   f.open( QFile::ReadOnly );
   QJsonDocument doc = QJsonDocument::fromJson( f.readAll() );

   GlobalSymmetry::setSymmetry( doc["symmetry"].toString().toStdString() );
   MatrixIndexMap::update();

   for ( const QJsonValue& vertex_ : doc["vertices"].toArray() )
   {
      QJsonObject vertex = vertex_.toObject();
      dual->addVertex( vertex["color"].toInt(), XYZ( vertex["x"].toDouble(), vertex["y"].toDouble(), vertex["z"].toDouble() ) );      
      dual->_Vertices.back()._SymmetryMap = MatrixSymmetryMap::symmetryFor( dual->_Vertices.back()._Pos );
   }
   for ( const QJsonValue& edge_ : doc["edges"].toArray() )
   {
      QJsonArray edge = edge_.toArray();
      dual->toggleEdge( Dual::VertexPtr( edge[0].toInt(), QMtx4x4() ), Dual::VertexPtr( edge[1].toInt(), MatrixIndexMap::at( edge[2].toInt() ) ), true/*only add edges*/ );
   }
   return dual;
}



SphereColoring::SphereColoring( QWidget *parent )
   : QMainWindow( parent )
{
   ui.setupUi( this );

   double radius = 5.0;
   //shared_ptr<Dual> dual = makeDual( radius );
   //shared_ptr<Dual> dual = loadDual( "GP7_1.dual" );
   //shared_ptr<Dual> dual = loadDual( "empty.dual" );
   shared_ptr<Dual> dual = loadDual( "temp.dual" );
   //shared_ptr<Dual> dual = loadDual( "GP10_6!.dual" );
   radius = dual->_Vertices[0]._Pos.len();
   //shared_ptr<Graph> graph = makeGraph( dual, radius );
   shared_ptr<Graph> graph = nullptr;
   _Simulation.init( dual, graph, radius );
   ui.drawing->_Simulation = &_Simulation;
   ui.lineEdit0->setText( "10" );
   ui.lineEdit1->setText( QString::number( _Simulation._Radius ) );

   connect( ui.yRotationSlider, &QSlider::valueChanged, [this]( int ) {
      ui.drawing->_YRotation = ui.yRotationSlider->value() * 1.;
      redrawSim();
   } );
   connect( ui.xRotationSlider, &QSlider::valueChanged, [this]( int ) {
      ui.drawing->_XRotation = ui.xRotationSlider->value() * 1.;
      redrawSim();
   } );
   connect( ui.zoomSlider, &QSlider::valueChanged, [this]( int ) {
      ui.drawing->_Zoom = 1. + ui.zoomSlider->value()/100.*3.;
      redrawSim();
   } );

   connect( ui.lineEdit0, &QLineEdit::editingFinished, [this]() {
      ui.drawing->_Custom[0] = ui.lineEdit0->text().toDouble();
      redrawSim();
   } );
   connect( ui.lineEdit1, &QLineEdit::editingFinished, [this]() {
      //ui.drawing->_Custom[1] = ui.lineEdit1->text().toDouble();
      _Simulation._Radius = ui.lineEdit1->text().toDouble();  
      if ( _Simulation._Radius < .5 )
      {
         _Simulation._Radius = 5;
         ui.lineEdit1->setText( QString::number( _Simulation._Radius ) );
      }
      _Simulation.normalizeVertices();
      redrawSim();
   } );
   //connect( ui.stepButton, &QPushButton::pressed, [this]() {
   //   _Simulation.step();
   //   redrawSim();
   //} );
   connect( ui.playButton, &QPushButton::pressed, [this]() {
      if ( _Timer.isActive() )
      {
         _Timer.stop();
         ui.playButton->setText( "Play" );
      }
      else
      {
         _Timer.start( 50 );
         ui.playButton->setText( "Pause" );
      }
   } );
   connect( &_Timer, &QTimer::timeout, [this]() {
      double error = _Simulation.step( 500 );
      ui.errorLabel->setText( "Err:" + QString::number( error ) );
      ui.paddingErrorLabel->setText( "Pad:" + QString::number( _Simulation._PaddingError ) );      
      redrawSim();
   } );
   
   connect( ui.showDualCheckBox      , &QCheckBox::toggled, [this]() { ui.drawing->_ShowDual       = ui.showDualCheckBox      ->isChecked(); redrawSim(); } );
   connect( ui.drawCurvesCheckBox    , &QCheckBox::toggled, [this]() { ui.drawing->_DrawCurves     = ui.drawCurvesCheckBox    ->isChecked(); redrawSim(); } );
   connect( ui.drawRigidsCheckBox    , &QCheckBox::toggled, [this]() { ui.drawing->_DrawRigidEDs   = ui.drawRigidsCheckBox    ->isChecked(); redrawSim(); } );
   connect( ui.labelVerticesCheckBox , &QCheckBox::toggled, [this]() { ui.drawing->_LabelVertices  = ui.labelVerticesCheckBox ->isChecked(); redrawSim(); } );
   connect( ui.labelTilesCheckBox    , &QCheckBox::toggled, [this]() { ui.drawing->_LabelTiles     = ui.labelTilesCheckBox    ->isChecked(); redrawSim(); } );
   connect( ui.drawSectors           , &QCheckBox::toggled, [this]() { ui.drawing->_DrawSectors    = ui.drawSectors           ->isChecked(); redrawSim(); } );
   connect( ui.showViolationsCheckBox, &QCheckBox::toggled, [this]() { ui.drawing->_ShowViolations = ui.showViolationsCheckBox->isChecked(); redrawSim(); } );

   connect( ui.dualToGraphButton, &QPushButton::pressed, [this]() {
      shared_ptr<Graph> graph = makeGraph( _Simulation._Dual, _Simulation._Radius );
      _Simulation.init( _Simulation._Dual, graph, _Simulation._Radius );
      redrawSim();
   } );


   connect( ui.loadButton, &QPushButton::pressed, [this]() {
      QString filename = QFileDialog::getOpenFileName( this, "Load Graph", QString(), "*.dual" );
      shared_ptr<Dual> dual = loadDual( filename );
      if ( !dual )
         return;
      _Simulation._Dual = dual;
      _Simulation._Radius = dual->_Vertices[0]._Pos.len();
      ui.lineEdit1->setText( QString::number( _Simulation._Radius ) );

      ui.dualToGraphButton->click();
      //redrawSim();
   } );
   connect( ui.saveButton, &QPushButton::pressed, [this]() {
      QString filename = QFileDialog::getSaveFileName( this, "Save Graph", QString(), "*.dual" );
      saveDual( filename, *_Simulation._Dual );

      //for ( const Graph::Vertex& v : _Simulation._Graph->_Vertices )
      //{
      //   char buf[500];
      //   sprintf( buf, "%2d: %15.12f %15.12f %15.12f", v._Index, v._Pos.x, v._Pos.y, v._Pos.z );
      //   qDebug().noquote() << buf;
      //}
      //for ( auto config : _Simulation._Graph->_Ico.matrices() )
      //{
      //   char buf[500];
      //   qDebug() << "";
      //   sprintf( buf, "%15.12f %15.12f %15.12f", config.m[0].x, config.m[1].x, config.m[2].x );
      //   qDebug().noquote() << buf;
      //   sprintf( buf, "%15.12f %15.12f %15.12f", config.m[0].y, config.m[1].y, config.m[2].y );
      //   qDebug().noquote() << buf;
      //   sprintf( buf, "%15.12f %15.12f %15.12f", config.m[0].z, config.m[1].z, config.m[2].z );
      //   qDebug().noquote() << buf;
      //}
      //vector<QMtx4x4> v = {
      //_Simulation._Graph->_Ico.map( {0,1,2}, {3,0,2} ),
      //_Simulation._Graph->_Ico.map( {0,1,2}, {7,6,8} ),
      //_Simulation._Graph->_Ico.map( {0,1,2}, {1,0,5} ),
      //_Simulation._Graph->_Ico.map( {0,1,2}, {1,2,0} ), };
      //for ( auto m : v )
      //{
      //   char buf[500];
      //   qDebug() << "";
      //   sprintf( buf, "%15.12f %15.12f %15.12f", m[0].x, m[1].x, m[2].x );
      //   qDebug().noquote() << buf;
      //   sprintf( buf, "%15.12f %15.12f %15.12f", m[0].y, m[1].y, m[2].y );
      //   qDebug().noquote() << buf;
      //   sprintf( buf, "%15.12f %15.12f %15.12f", m[0].z, m[1].z, m[2].z );
      //   qDebug().noquote() << buf;
      //}
   } );


   connect( ui.drawing, &Drawing::press, [this]( QMouseEvent* event ) {
      if ( event->buttons().testFlag( Qt::LeftButton ) )
         handleMouse( event->pos(), true, true, false );
      if ( event->buttons().testFlag( Qt::RightButton ) )
         handleRightButton( event->pos(), true, true, false );
   } );

   connect( ui.drawing, &Drawing::move, [this]( QMouseEvent* event ) {
      handleMouse( event->pos(), true, false, false );
      if ( event->buttons().testFlag( Qt::RightButton ) )
         handleRightButton( event->pos(), true, false, false );
   } );

   connect( ui.drawing, &Drawing::release, [this]( QMouseEvent* event ) {
      //if ( event->buttons().testFlag( Qt::LeftButton ) )
         handleMouse( event->pos(), false, false, true );
      //if ( event->buttons().testFlag( Qt::RightButton ) )
         handleRightButton( event->pos(), false, false, true );
   } );
   

   QObject::connect( new QShortcut(QKeySequence(Qt::Key_0), this ), &QShortcut::activated, [this]() { addVertex( 0 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_R), this ), &QShortcut::activated, [this]() { addVertex( 0 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_1), this ), &QShortcut::activated, [this]() { addVertex( 1 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_B), this ), &QShortcut::activated, [this]() { addVertex( 1 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_2), this ), &QShortcut::activated, [this]() { addVertex( 2 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_Y), this ), &QShortcut::activated, [this]() { addVertex( 2 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_3), this ), &QShortcut::activated, [this]() { addVertex( 3 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_P), this ), &QShortcut::activated, [this]() { addVertex( 3 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_4), this ), &QShortcut::activated, [this]() { addVertex( 4 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_G), this ), &QShortcut::activated, [this]() { addVertex( 4 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_5), this ), &QShortcut::activated, [this]() { addVertex( 5 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_O), this ), &QShortcut::activated, [this]() { addVertex( 5 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_6), this ), &QShortcut::activated, [this]() { addVertex( 6 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_C), this ), &QShortcut::activated, [this]() { addVertex( 6 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_7), this ), &QShortcut::activated, [this]() { addVertex( 7 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_W), this ), &QShortcut::activated, [this]() { addVertex( 7 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_8), this ), &QShortcut::activated, [this]() { addVertex( 8 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_K), this ), &QShortcut::activated, [this]() { addVertex( 8 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_9), this ), &QShortcut::activated, [this]() { addVertex( BLANK_COLOR ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_N), this ), &QShortcut::activated, [this]() { addVertex( BLANK_COLOR ); } );

   QObject::connect( new QShortcut(QKeySequence(Qt::Key_F1), this ), &QShortcut::activated, [this]() { toggleSymmetryVertex( 0 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_F2), this ), &QShortcut::activated, [this]() { toggleSymmetryVertex( 1 ); } );
   QObject::connect( new QShortcut(QKeySequence(Qt::Key_F3), this ), &QShortcut::activated, [this]() { toggleSymmetryVertex( 2 ); } );

   QObject::connect( new QShortcut(QKeySequence(Qt::Key_Delete), this ), &QShortcut::activated, [this]() { deleteVertex(); } );


   //connect( ui.permSlider, &QSlider::valueChanged, [this]( int ) {
   //   ui.drawing->_PermIndex = ui.permSlider->value();
   //   redrawSim();
   //} );
   redrawSim();
}

void SphereColoring::handleRightButton( const QPoint& mousePos, bool isMove, bool isClick, bool isUnclick )
{
   if ( isClick )
   {
      _MouseRightButtonDownPos = mousePos;
      _PreDragModelRotation = ui.drawing->_ModelRotation;
   }
   if ( isMove )
   {
      QPointF p = mousePos - _MouseRightButtonDownPos;
      XYZ axis = p == QPointF(0,0) ? XYZ(0,0,1) : XYZ( -p.y(), -p.x(), 0 );
      ui.drawing->_ModelRotation = QMtx4x4::rotation( axis, QLineF( QPointF(), p ).length() * 1.99 / ui.drawing->_Zoom / ui.drawing->height() ) * _PreDragModelRotation;
      redrawSim();
   }
}

void SphereColoring::redrawSim() const 
{ 
   ui.drawing->_Radius = _Simulation._Radius;
   ui.drawing->updateLabel(); 
}

void SphereColoring::addVertex( int color )
{
   //if ( color >= 7 )
   //   return;

   XYZ modelPos;
   
   Dual::VertexPtr existingVertex = ui.drawing->dualVertexNearest( ui.drawing->mousePos(), 6. );
   if ( existingVertex.isValid() )
   {
      _Simulation._Dual->setColorOf( existingVertex, color );
   }
   else
   {
      if ( ui.drawing->bitmapToModel( ui.drawing->mousePos(), modelPos ) );
         _Simulation._Dual->addVertex( color, modelPos );
   }

   redrawSim();
}

void SphereColoring::handleMouse( const QPoint& mousePos, bool isMove, bool isClick, bool isUnclick )
{
   if ( isUnclick )
   {
      if ( isKeyDown( 'E' ) )
      {
         _Simulation._Dual->toggleEdge( _EdgeSelectVtx, ui.drawing->dualVertexNearest( mousePos, 30. ) );
         redrawSim();
      }

      _DragDualVtx = Dual::VertexPtr();
      _EdgeSelectVtx = Dual::VertexPtr();
   }

   if ( isKeyDown( 'M' ) && isClick )
   {
      _DragDualVtx = ui.drawing->dualVertexNearest( mousePos, 8. );
   }

   if ( isKeyDown( 'E' ) && isClick )
   {
      _EdgeSelectVtx = ui.drawing->dualVertexNearest( mousePos, 12. );
   }

   if ( isKeyDown( 'M' ) && isMove )
   {
      XYZ modelPos;
      if ( _DragDualVtx.isValid() && ui.drawing->bitmapToModel( mousePos, modelPos ) )
         _Simulation._Dual->setPos( _DragDualVtx, modelPos );
      redrawSim();
   }   

   if ( isKeyDown( 'E' ) && isMove )
   {
      _DragDualVtx = ui.drawing->dualVertexNearest( mousePos, 50. );
   }
}

void SphereColoring::toggleSymmetryVertex( int idx )
{
   vector<XYZ> sym = GlobalSymmetry::symmetryPoints( _Simulation._Radius );
   if ( idx >= (int) sym.size() )
      return;

   XYZ p = sym[idx].normalized() * _Simulation._Radius;

   int vertexIdx = -1;
   for ( const Dual::Vertex& vtx : _Simulation._Dual->_Vertices )
   {
      if ( vtx._Pos.dist2( p ) < 1e-12 )
         vertexIdx = vtx._Index;
   }

   if ( vertexIdx >= 0 )
   {
      _Simulation._Dual->deleteVertex( vertexIdx );
   }
   else
   {
      _Simulation._Dual->addVertex( 0, p );
      _Simulation._Dual->_Vertices.back()._SymmetryMap = MatrixSymmetryMap::symmetryFor( p );
   }

   redrawSim();
}


void SphereColoring::deleteVertex()
{
   Dual::VertexPtr existingVertex = ui.drawing->dualVertexNearest( ui.drawing->mousePos(), 6. );
   if ( existingVertex.isValid() )
   {
      _Simulation._Dual->deleteVertex( existingVertex._Index );
   }

   redrawSim();
}
