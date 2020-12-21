#pragma once

#include "Model.h"
#include <memory>

class Simulation
{
public:
   void init( shared_ptr<Dual> dual, std::shared_ptr<Graph> graph, double radius );
   void normalizeVertices();
   double step( double& paddingError );
   double step( int numSteps );

public:
   double _Radius = 1;
   double _Padding = .0001;
   double _PaddingError = 0;
   shared_ptr<Graph> _Graph;
   vector<Graph::KeepCloseFar> _KeepCloseFars;
   vector<Graph::LineVertexConstraint> _LineVertexConstraints;
   shared_ptr<Dual> _Dual;
};

