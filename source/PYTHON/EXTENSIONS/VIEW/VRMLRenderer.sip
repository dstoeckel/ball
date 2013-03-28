class VRMLRenderer : SceneExporter
{
%TypeHeaderCode
	#include <BALL/VIEW/RENDERING/RENDERERS/VRMLRenderer.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  public:

	VRMLRenderer(const String&);
	~VRMLRenderer();
	virtual void clear();
	String VRMLColorRGBA(const ColorRGBA&);
	String VRMLVector3(Vector3);
	void VRMLColor(const ColorRGBA&);
	virtual bool init(const Stage*, float, float);
	bool exportOneRepresentation(const Representation*);
};
