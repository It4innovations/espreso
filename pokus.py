
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.1.0-136-g1cedb6a 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.1.0-136-g1cedb6a

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1192, 600]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0, 0, 0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0, 0, 0]
      renderView1.CameraFocalPoint = [14.999999821186066, 14.999999821186066, 14.999999821186066]
      renderView1.CameraViewUp = [0, 0, 0]
      renderView1.CameraParallelScale = 25
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=1, magnification=1, width=1192, height=600, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      results0vtk = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Warp By Vector'
      warpByVector1 = WarpByVector(Input=results0vtk)
      warpByVector1.Vectors = ['POINTS', 'displacement']
      warpByVector1.ScaleFactor = 16390.0

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'displacement'
      displacementLUT = GetColorTransferFunction('displacement')
      displacementLUT.RGBPoints = [5.868136485274054e-16, 0.231373, 0.298039, 0.752941, 8.099107620490154e-05, 0.865003, 0.865003, 0.865003, 0.0001619821524092163, 0.705882, 0.0156863, 0.14902]
      displacementLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'displacement'
      displacementPWF = GetOpacityTransferFunction('displacement')
      displacementPWF.Points = [5.868136485274054e-16, 0.0, 0.5, 0.0, 0.0001619821524092163, 1.0, 0.5, 0.0]
      displacementPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from warpByVector1
      warpByVector1Display = Show(warpByVector1, renderView1)
      # trace defaults for the display properties.
      warpByVector1Display.ColorArrayName = ['POINTS', 'displacement']
      warpByVector1Display.LookupTable = displacementLUT
      warpByVector1Display.OSPRayScaleArray = 'displacementIsSet'
      warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      warpByVector1Display.GlyphType = 'Arrow'
      warpByVector1Display.ScalarOpacityUnitDistance = 2.265150654738395

      # show color legend
      warpByVector1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for displacementLUT in view renderView1
      displacementLUTColorBar = GetScalarBar(displacementLUT, renderView1)
      displacementLUTColorBar.Title = 'displacement'
      displacementLUTColorBar.ComponentTitle = 'Magnitude'

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(warpByVector1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1]}
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
coprocessor.EnableLiveVisualization(True, 1)


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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
