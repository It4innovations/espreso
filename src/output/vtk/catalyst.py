
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.0.1 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.0.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1790, 1196]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.5002305023372173, 0.5003451500087976, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.5002305023372173, 0.5003451500087976, 10000.0]
      renderView1.CameraFocalPoint = [0.5002305023372173, 0.5003451500087976, 0.0]
      renderView1.CameraParallelScale = 0.6311503291488461
      renderView1.Background = [0.32, 0.34, 0.43]

      # scale = datadescription.GetUserData().GetAbstractArray("scale").GetValue(0);
      # label = datadescription.GetUserData().GetAbstractArray("label").GetValue(0);
      # print scale
      # print label

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=1, magnification=1, width=1790, height=1196, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      results0vtk = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'temperatures0i0'
      temperatures0i0LUT = GetColorTransferFunction('temperatures0i0')
      temperatures0i0LUT.RGBPoints = [-4.77577e-06, 0.231373, 0.298039, 0.752941, 0.9999976121149999, 0.865003, 0.865003, 0.865003, 2.0, 0.705882, 0.0156863, 0.14902]
      temperatures0i0LUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'temperatures0i0'
      temperatures0i0PWF = GetOpacityTransferFunction('temperatures0i0')
      temperatures0i0PWF.Points = [-4.77577e-06, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
      temperatures0i0PWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from results0vtk
      results0vtkDisplay = Show(results0vtk, renderView1)
      # trace defaults for the display properties.
      results0vtkDisplay.ColorArrayName = ['POINTS', 'results']
      results0vtkDisplay.LookupTable = temperatures0i0LUT
      results0vtkDisplay.GlyphType = 'Arrow'
      # results0vtkDisplay.ScalarOpacityUnitDistance = 0.036909895158202304
      # results0vtkDisplay.SetScaleArray = ['POINTS', 'results']
      # results0vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      # results0vtkDisplay.OpacityArray = ['POINTS', 'results']
      # results0vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      results0vtkDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for temperatures0i0LUT in view renderView1
      temperatures0i0LUTColorBar = GetScalarBar(temperatures0i0LUT, renderView1)
      temperatures0i0LUTColorBar.Title = 'results'
      temperatures0i0LUTColorBar.ComponentTitle = ''
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=True)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
