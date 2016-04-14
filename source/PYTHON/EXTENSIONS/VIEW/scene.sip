// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id$
//

class Scene
  : ModularWidget
{
%TypeHeaderCode
	#include <BALL/VIEW/WIDGETS/scene.h>
	using namespace BALL;
	using namespace BALL::VIEW;
%End
  private:
    // Scene is not copyable
    Scene(const Scene&);
  public:
//  void dump(std::ostream& s = std::cout, Size depth = 0) const;
//  BALL_EMBEDDABLE(Scene,);

// methods implemented by the BALL_EMBEDDED macro
static Size countInstances();
static Scene* getInstance(Position);

  Scene();
//  Scene(QWidget*, const char*, WFlags);
//  const Scene& operator =  (const Scene&);
//  void fetchPreferences(INIFile&);
//  void finalizePreferencesTab(Preferences &preferences);
//  void finalizeWidget(MainControl&);
//  void initializePreferencesTab(Preferences &preferences);
//  void initializeWidget(MainControl&);
//  void set(const Scene&);
  ~Scene();

//////////////////	void addStructure(String name);
	void setCursor(String c);
	void setElementCursor(int);

	virtual void onNotify(Message*);
	virtual bool exportScene(SceneExporter&) const;
	virtual void initializeWidget(MainControl&);
	virtual void fetchPreferences(INIFile&);
	virtual void writePreferences(INIFile&);
	virtual void checkMenu(MainControl&);
	virtual bool isValid() const;
//	virtual void dump(std::ostream& = std::cout, Size = 0) const;
	Stage* getStage();
//	void setInfo(const QString& );
	Vector2 map3DToViewport(const Vector3& );
	Vector3 mapViewportTo3D(const Vector2& );
	Vector3 mapViewportTo3D(int , int );
//	void pickObjects(const QPoint& , list<GeometricObject*>& );
//	void pickObjects(const QPoint& , const QPoint& , list<GeometricObject*>& );
	void merge(Composite* , Composite* );
//	std::list<AtomContainer*> getContainers();
//	void drawLine(const QPointF& , const QPointF& , QPainter* );
//	void drawLine(const Vector2& , const Vector2& , QPainter* );
	void deselect(bool );
//	void getClickedItems(const QPoint& );

//	const Stage* getStage() const;
	void setCamera(const Camera&);
//	List<Camera>& getAnimationPoints();
	static void setAnimationSmoothness(float);
	static float getAnimationSmoothness();
	void setDefaultLighting(bool /*update_GL*/ = true);
	void initializePreferencesTab(Preferences&);
	void finalizePreferencesTab(Preferences&);
	virtual void applyPreferences();

//	static void setMouseSensitivity(float sensitivity);
//	static float getMouseSensitivity();
//	static void setMouseWheelSensitivity(float sensitivity);
//	static float getMouseWheelSensitivity();
	
//	static void setShowLightSources(bool state);
//	static bool showLightSourcesEnabled();
//	GLRenderer& getGLRenderer();
//	static bool stereoBufferSupportTest();
//	float getMousePositionX();
//	float getMousePositionY();
	bool exportPNG(const String&);
	virtual void setWidgetVisible(bool);
	void setOffScreenRendering(bool, Size);
//	virtual void keyPressEvent(QKeyEvent*);
//	virtual void mouseDoubleClickEvent(QMouseEvent*);
//	ModeType getMode() const;
//	virtual void setMode(ModeType);
//	static void setScreenShotNumber(Position pos);
//	static void setPOVNumber(Position pos);
	void rotate(float /*degree_right*/, float /*degree_up*/);
	void rotateClockwise(float /*degree*/);
	void move(Vector3);
	void moveComposites(PyCompositeList, Vector3);
	void rotateComposites(PyCompositeList, float /*degree_right*/, float /*degree_up*/, float de/*gree_clockwise*/ = 0);
	bool isAnimationRunning() const;
	void setTurnPoint(const Vector3&);
	const Vector3& getTurnPoint() const;
	void setFullScreen(bool);
	void setFPSEnabled(bool);
	void showText(const String&, Size /*font_size*/ = 20);
//	virtual void addToolBarEntries(QToolBar*);

	public slots:

	void createCoordinateSystem();
	void createCoordinateSystemAtOrigin();
	String exportPNG();
	void exportPOVRay();
	void exportNextPOVRay();
	void printScene();
	void showExportPNGDialog();
	void showExportVRMLDialog();
	virtual void switchShowWidget();
	void exitStereo();
	void enterStereo();
	void enterActiveStereo();
	void enterDualStereo();
	void clearRecordedAnimation();
	void startAnimation();
	void stopAnimation();
	void switchToLastMode();
	void switchShowGrid();
	void showInfos();
	void setupViewVolume();
	void createNewMolecule();
	void saturateWithHydrogens();
	void optimizeStructure();
	void setEditElementType(int element_number);
	int getEditElementType();
	void toggleContinuousLoop();
	void startContinuousLoop();
	void stopContinuousLoop();
	void applyStereoDefaults();
};
