
PROJECT( SmartVesselExtraction )

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )


SET(ITK_DIR C:/VITK_10/ITK3201/Build)
SET(VTK_DIR C:/VITK_10/VTK580/Build)




# Find TIK
FIND_PACKAGE( ITK REQUIRED )
IF( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ENDIF( ITK_FOUND )

FIND_PACKAGE(VTK REQUIRED)
IF(NOT VTK_USE_RENDERING)
  MESSAGE(FATAL_ERROR "Example ${PROJECT_NAME} requires VTK_USE_RENDERING.")
ENDIF(NOT VTK_USE_RENDERING)
INCLUDE(${VTK_USE_FILE})


ADD_EXECUTABLE( SmartVesselExtraction 
MainSmartVesselExtraction.cxx
TestVesselExtraction.cxx

)

# Link to some libraries
TARGET_LINK_LIBRARIES( SmartVesselExtraction 
${ITK_LIBRARIES}

MapReduceMPI.lib
mpistubs.lib
vtkalglib.lib
vtkCharts.lib
vtkCommon.lib
vtkDICOMParser.lib
vtkexoIIc.lib
vtkexpat.lib
vtkFiltering.lib
vtkfreetype.lib
vtkftgl.lib
vtkGenericFiltering.lib
vtkGeovis.lib
vtkGraphics.lib
vtkhdf5.lib
vtkHybrid.lib
vtkImaging.lib
vtkInfovis.lib
vtkIO.lib
vtkjpeg.lib
vtklibxml2.lib
vtkmetaio.lib
vtkNetCDF.lib
vtkNetCDF_cxx.lib
vtkpng.lib
vtkproj4.lib
vtkRendering.lib
vtksqlite.lib
vtksys.lib
vtktiff.lib
vtkverdict.lib
vtkViews.lib
vtkVolumeRendering.lib
vtkWidgets.lib
vtkzlib.lib
)
