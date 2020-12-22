#pragma once

#include <QWidget>
#include "ui_Drawing.h"
#include <memory>
#include <set>
#include "Model.h"

class Graph;
class Simulation;

class Drawing : public QWidget
{
   Q_OBJECT

public:
   Drawing( QWidget *parent = Q_NULLPTR );
   ~Drawing();

   QImage makeImage( Graph& graph );
   void updateLabel();

   QPoint mousePos() const;
   QMtx4x4 modelToBitmap() const;
   QMtx4x4 modelToBitmapNoRot() const;
   QMtx4x4 modelRotation() const;

   bool bitmapToModel( const QPoint& p, XYZ& modelPos ) const;
   bool isOnNearSide( const XYZ& p ) const;
   Dual::VertexPtr dualVertexNearest( const QPoint& mousePos, double maxPixelDist );

private:
   void resizeEvent( QResizeEvent *event ) override;
   //void mousePressEvent(QMouseEvent * event);
   
   void mousePressEvent( QMouseEvent * event ) override { emit press( event ); debugClick( event ); }
   void mouseReleaseEvent( QMouseEvent * event ) override { emit release( event ); }
   void mouseMoveEvent( QMouseEvent * event ) override { emit move( event ); }

   void debugClick( QMouseEvent* event );

signals:
   void press( QMouseEvent * event );
   void release( QMouseEvent * event );
   void move( QMouseEvent * event );

private:
   Ui::Drawing ui;

public:
   QMtx4x4 _ModelRotation;
   double _Radius = 0;
   double _YRotation = 0;
   double _XRotation = 0;
   double _Zoom = 1.;
   double _Custom[10] = { 0 };
   QPoint _ClickPos;
   bool _ShowDual = true;
   bool _DrawCurves = true;
   bool _DrawRigidEDs = true;
   bool _LabelVertices = false;
   bool _LabelTiles = false;
   bool _DrawSectors = true;
   bool _ShowViolations = true;

   bool _DrawZoneOfExclusions = true;
   const Simulation* _Simulation;
   std::vector<Graph::VertexPtr> _SelectedVertices;
};

