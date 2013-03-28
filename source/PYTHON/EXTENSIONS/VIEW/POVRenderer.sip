// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id$
//

class POVRenderer : SceneExporter
{
%TypeHeaderCode
	#include <BALL/VIEW/RENDERING/RENDERERS/POVRenderer.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:
  POVRenderer(const String&);
//  String POVColorRGBA(const ColorRGBA&);
//  String POVFinish(const String&, const ColorRGBA&);
//  String POVVector3(Vector3);
  void clear();

	bool exportOneRepresentation(const Representation*);
	virtual bool init(const Stage*, float, float);
  ~POVRenderer();
};
