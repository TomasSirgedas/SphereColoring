#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_SphereColoring.h"

#include "Simulation.h"
#include <QTimer>

class SphereColoring : public QMainWindow
{
   Q_OBJECT

public:
   SphereColoring( QWidget *parent = Q_NULLPTR );

   void redrawSim() const;
   void addVertex( int color );   

   void handleMouse( const QPoint& mousePos, bool isMove, bool isClick, bool isUnclick );
   void handleRightButton( const QPoint& mousePos, bool isMove, bool isClick, bool isUnclick );
   
   Dual::VertexPtr dualVertexNearest( const QPoint& mousePos );

private:
   Ui::SphereColoringClass ui;

   Simulation _Simulation;
   QTimer _Timer;
   Dual::VertexPtr _DragDualVtx;
   Dual::VertexPtr _EdgeSelectVtx;

   QPoint _MouseRightButtonDownPos;
   QMtx4x4 _PreDragModelRotation;
};
