### list all filenames of the directory here ###
SET(SOURCES_LIST
	renderer.C
	POVRenderer.C
	STLRenderer.C
	tilingRenderer.C
	VRMLRenderer.C
	XML3DRenderer.C
)

ADD_VIEW_SOURCES("VIEW/RENDERING/RENDERERS" "${SOURCES_LIST}")
