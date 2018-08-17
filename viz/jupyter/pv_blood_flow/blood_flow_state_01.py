# state file generated using paraview version 5.5.2

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

import sys


start_frame=int(sys.argv[1])
num_frames=int(sys.argv[2])
DEST_IMAGE_DIR=sys.argv[3]

DATA_DIR="/projects/ATPESC2018/VIS_DATA/BLOODFLOW_ANIMATION_DATA"

RBC_GOOD_DATA_FILES = []
RBC_BAD_DATA_FILES = []

for i in range(100):
    temp_name = "%s/rbc_good_step_%04d.vtp" % (DATA_DIR, i)
    RBC_GOOD_DATA_FILES.append(temp_name)
    temp_name = "%s/rbc_bad_step_%04d.vtp" % (DATA_DIR, i)
    RBC_BAD_DATA_FILES.append(temp_name)


# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [960, 540]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [19.630000710487366, -0.3942399024963379, -0.927609920501709]
renderView1.StereoType = 0
renderView1.CameraPosition = [32.993379066782744, -20.602286555938573, -52.17305711591688]
renderView1.CameraFocalPoint = [20.551871526660904, -1.7882892188435457, -4.462770040319944]
renderView1.CameraViewUp = [0.18596694126457572, -0.8965400471506745, 0.40203512360459787]
renderView1.CameraParallelScale = 24.19731677878087
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML PolyData Reader'
rbc_good_step_0000vtp = XMLPolyDataReader(FileName=RBC_GOOD_DATA_FILES)
rbc_good_step_0000vtp.PointArrayStatus = ['Velocity', 'Normals']

# create a new 'XML PolyData Reader'
rbc_bad_step_0000vtp = XMLPolyDataReader(FileName=RBC_BAD_DATA_FILES)
rbc_bad_step_0000vtp.PointArrayStatus = ['Velocity', 'Normals']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from rbc_bad_step_0000vtp
rbc_bad_step_0000vtpDisplay = Show(rbc_bad_step_0000vtp, renderView1)

# trace defaults for the display properties.
rbc_bad_step_0000vtpDisplay.Representation = 'Surface'
rbc_bad_step_0000vtpDisplay.ColorArrayName = [None, '']
rbc_bad_step_0000vtpDisplay.DiffuseColor = [0.7686274509803922, 0.0, 0.0]
rbc_bad_step_0000vtpDisplay.Specular = 0.46
rbc_bad_step_0000vtpDisplay.SpecularPower = 54.0
rbc_bad_step_0000vtpDisplay.Ambient = 0.14
rbc_bad_step_0000vtpDisplay.OSPRayScaleArray = 'Normals'
rbc_bad_step_0000vtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
rbc_bad_step_0000vtpDisplay.SelectOrientationVectors = 'Velocity'
rbc_bad_step_0000vtpDisplay.ScaleFactor = 4.5540001630783085
rbc_bad_step_0000vtpDisplay.SelectScaleArray = 'None'
rbc_bad_step_0000vtpDisplay.GlyphType = 'Arrow'
rbc_bad_step_0000vtpDisplay.GlyphTableIndexArray = 'None'
rbc_bad_step_0000vtpDisplay.GaussianRadius = 0.22770000815391542
rbc_bad_step_0000vtpDisplay.SetScaleArray = ['POINTS', 'Normals']
rbc_bad_step_0000vtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
rbc_bad_step_0000vtpDisplay.OpacityArray = ['POINTS', 'Normals']
rbc_bad_step_0000vtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
rbc_bad_step_0000vtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
rbc_bad_step_0000vtpDisplay.SelectionCellLabelFontFile = ''
rbc_bad_step_0000vtpDisplay.SelectionPointLabelFontFile = ''
rbc_bad_step_0000vtpDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
rbc_bad_step_0000vtpDisplay.ScaleTransferFunction.Points = [-0.9994208812713623, 0.0, 0.5, 0.0, 0.9991291165351868, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
rbc_bad_step_0000vtpDisplay.OpacityTransferFunction.Points = [-0.9994208812713623, 0.0, 0.5, 0.0, 0.9991291165351868, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
rbc_bad_step_0000vtpDisplay.DataAxesGrid.XTitleFontFile = ''
rbc_bad_step_0000vtpDisplay.DataAxesGrid.YTitleFontFile = ''
rbc_bad_step_0000vtpDisplay.DataAxesGrid.ZTitleFontFile = ''
rbc_bad_step_0000vtpDisplay.DataAxesGrid.XLabelFontFile = ''
rbc_bad_step_0000vtpDisplay.DataAxesGrid.YLabelFontFile = ''
rbc_bad_step_0000vtpDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
rbc_bad_step_0000vtpDisplay.PolarAxes.PolarAxisTitleFontFile = ''
rbc_bad_step_0000vtpDisplay.PolarAxes.PolarAxisLabelFontFile = ''
rbc_bad_step_0000vtpDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
rbc_bad_step_0000vtpDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from rbc_good_step_0000vtp
rbc_good_step_0000vtpDisplay = Show(rbc_good_step_0000vtp, renderView1)

# trace defaults for the display properties.
rbc_good_step_0000vtpDisplay.Representation = 'Surface'
rbc_good_step_0000vtpDisplay.ColorArrayName = [None, '']
rbc_good_step_0000vtpDisplay.DiffuseColor = [0.0, 0.6666666666666666, 1.0]
rbc_good_step_0000vtpDisplay.Specular = 0.34
rbc_good_step_0000vtpDisplay.SpecularPower = 83.0
rbc_good_step_0000vtpDisplay.Ambient = 0.17
rbc_good_step_0000vtpDisplay.OSPRayScaleArray = 'Normals'
rbc_good_step_0000vtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
rbc_good_step_0000vtpDisplay.SelectOrientationVectors = 'Velocity'
rbc_good_step_0000vtpDisplay.ScaleFactor = 4.37199981212616
rbc_good_step_0000vtpDisplay.SelectScaleArray = 'None'
rbc_good_step_0000vtpDisplay.GlyphType = 'Arrow'
rbc_good_step_0000vtpDisplay.GlyphTableIndexArray = 'None'
rbc_good_step_0000vtpDisplay.GaussianRadius = 0.218599990606308
rbc_good_step_0000vtpDisplay.SetScaleArray = ['POINTS', 'Normals']
rbc_good_step_0000vtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
rbc_good_step_0000vtpDisplay.OpacityArray = ['POINTS', 'Normals']
rbc_good_step_0000vtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
rbc_good_step_0000vtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
rbc_good_step_0000vtpDisplay.SelectionCellLabelFontFile = ''
rbc_good_step_0000vtpDisplay.SelectionPointLabelFontFile = ''
rbc_good_step_0000vtpDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
rbc_good_step_0000vtpDisplay.ScaleTransferFunction.Points = [-0.9997668266296387, 0.0, 0.5, 0.0, 0.9994969964027405, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
rbc_good_step_0000vtpDisplay.OpacityTransferFunction.Points = [-0.9997668266296387, 0.0, 0.5, 0.0, 0.9994969964027405, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
rbc_good_step_0000vtpDisplay.DataAxesGrid.XTitleFontFile = ''
rbc_good_step_0000vtpDisplay.DataAxesGrid.YTitleFontFile = ''
rbc_good_step_0000vtpDisplay.DataAxesGrid.ZTitleFontFile = ''
rbc_good_step_0000vtpDisplay.DataAxesGrid.XLabelFontFile = ''
rbc_good_step_0000vtpDisplay.DataAxesGrid.YLabelFontFile = ''
rbc_good_step_0000vtpDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
rbc_good_step_0000vtpDisplay.PolarAxes.PolarAxisTitleFontFile = ''
rbc_good_step_0000vtpDisplay.PolarAxes.PolarAxisLabelFontFile = ''
rbc_good_step_0000vtpDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
rbc_good_step_0000vtpDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(rbc_bad_step_0000vtp)
# ----------------------------------------------------------------

tvals = rbc_good_step_0000vtp.TimestepValues

for i in range(start_frame, start_frame+num_frames):
    SetActiveView(renderView1)
    renderView1.ViewTime=tvals[i]

    IMAGE_FILENAME="%s/frame_%04d.png" % (DEST_IMAGE_DIR, i)
    print("Save: %s" % IMAGE_FILENAME)
    WriteImage(IMAGE_FILENAME)
