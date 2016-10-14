
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.3.1-805-g7670e63 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.3.1-805-g7670e63

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView2 = CreateView('RenderView')
      renderView2.ViewSize = [1268, 1001]
      renderView2.AxesGrid = 'GridAxes3DActor'
      renderView2.CenterOfRotation = [15.0, 15.0, 15.0]
      renderView2.StereoType = 0
      renderView2.CameraPosition = [98.52626368708493, -76.30116325299662, -30.50288786032776]
      renderView2.CameraFocalPoint = [16.662197115466938, 13.062278347419227, 21.00379719143666]
      renderView2.CameraViewUp = [-0.3333254352084389, 0.2226836782305292, -0.9161365256842602]
      renderView2.CameraParallelScale = 34.6254401549066
      renderView2.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView2,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1268, height=1001, cinema={})

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'PVTrivialProducer'
      # create a producer from a simulation input
      extractPVTrivialProducer1 = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Warp By Vector'
      warpByVector1 = WarpByVector(Input=extractPVTrivialProducer1)
      warpByVector1.Vectors = ['POINTS', 'displacement']
      warpByVector1.ScaleFactor = 1.0

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'displacement'
      displacementLUT = GetColorTransferFunction('displacement')
      # displacementLUT.RGBPoints = [5.028988747889835e-18, 0.231373, 0.298039, 0.752941, 7.472436195608509e-13, 0.865003, 0.865003, 0.865003, 1.494482210132954e-12, 0.705882, 0.0156863, 0.14902]
      displacementLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'displacement'
      displacementPWF = GetOpacityTransferFunction('displacement')
      displacementPWF.Points = [5.028988747889835e-18, 0.0, 0.5, 0.0, 1.494482210132954e-12, 1.0, 0.5, 0.0]
      displacementPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView2'
      # ----------------------------------------------------------------

      # show data from warpByVector1
      warpByVector1Display = Show(warpByVector1, renderView2)
      # trace defaults for the display properties.
      warpByVector1Display.ColorArrayName = ['POINTS', 'displacement']
      warpByVector1Display.LookupTable = displacementLUT
      warpByVector1Display.ScalarOpacityUnitDistance = 1.6193839689753242

      # show color legend
      warpByVector1Display.SetScalarBarVisibility(renderView2, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for displacementLUT in view renderView2
      displacementLUTColorBar = GetScalarBar(displacementLUT, renderView2)
      displacementLUTColorBar.Title = 'displacement'
      displacementLUTColorBar.ComponentTitle = 'Magnitude'
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
