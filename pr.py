
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
      renderView1.CenterOfRotation = [7.500000089406967, 7.500000089406967, 16.316250026226044]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [21.990637012885976, 69.82643606296968, 34.42881329060726]
      renderView1.CameraFocalPoint = [7.500000089406967, 7.500000089406967, 16.316250026226044]
      renderView1.CameraViewUp = [-0.1933027008293548, 0.31497678967077813, -0.9292059447833765]
      renderView1.CameraParallelScale = 17.212202322910283
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='aaaa.png', freq=1, fittoscreen=1, magnification=1, width=1192, height=600, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Unstructured Grid Reader'
      # create a producer from a simulation input
      results0vtu = coprocessor.CreateProducer(datadescription, 'input')
      
      # create a new 'Warp By Vector'
      scalef=datadescription.GetUserData().GetAbstractArray("Scale").GetValue(0);
      label=datadescription.GetUserData().GetAbstractArray("Label").GetValue(0);
      print(scalef,label)
      warpByVector1 = WarpByVector(Input=results0vtu)
      warpByVector1.Vectors = ['POINTS', label]
      warpByVector1.ScaleFactor = scalef

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'displacement'
      displacementLUT = GetColorTransferFunction(label)
      displacementLUT.RGBPoints = [1.0084477420735396e-16, 0.231373, 0.298039, 0.752941, 8.077062035316055e-05, 0.865003, 0.865003, 0.865003, 0.00016154124070622027, 0.705882, 0.0156863, 0.14902]
      displacementLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'displacement'
      displacementPWF = GetOpacityTransferFunction(label)
      displacementPWF.Points = [1.0084477420735396e-16, 0.0, 0.5, 0.0, 0.00016154124070622027, 1.0, 0.5, 0.0]
      displacementPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from warpByVector1
      warpByVector1Display = Show(warpByVector1, renderView1)
      # trace defaults for the display properties.
      warpByVector1Display.ColorArrayName = ['POINTS', label]
      warpByVector1Display.LookupTable = displacementLUT
      warpByVector1Display.OSPRayScaleArray = label
      warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      warpByVector1Display.GlyphType = 'Arrow'
      warpByVector1Display.ScalarOpacityUnitDistance = 3.442440464582057

      # show color legend
      warpByVector1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for displacementLUT in view renderView1
      displacementLUTColorBar = GetScalarBar(displacementLUT, renderView1)
      displacementLUTColorBar.Title = label
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
    #coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
