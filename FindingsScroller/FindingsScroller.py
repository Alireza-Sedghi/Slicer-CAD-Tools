import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

#
# FindingsScroller
#

class FindingsScroller(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Prostate Findings Scroller" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Alireza Mehrtash (BWH, UBC.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# FindingsScrollerWidget
#

class FindingsScrollerWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def onReload(self, moduleName="FindingsScroller"):
      """Generic reload method for any scripted module.
      ModuleWizard will subsitute correct default moduleName.
      """
      globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Reload and Test area
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)
    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "Freehand3DUltrasound Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)
    # uncomment the following line for debug/development.
    self.layout.addWidget(reloadCollapsibleButton)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    self.jsonFindingsPath = ctk.ctkPathLineEdit()
    self.jsonFindingsPath.currentPath = '/home/mehrtash/github/ProstateX/notebooks/findings.json'
    # self.dockerPath.setMaximumWidth(300)
    parametersFormLayout.addRow("Findings File Path:", self.jsonFindingsPath)

    #
    # Load Button
    #
    self.loadButton = qt.QPushButton("Load")
    self.loadButton.toolTip = "Run the algorithm."
    self.loadButton.enabled = False
    parametersFormLayout.addRow(self.loadButton)

    #
    # finding index 
    #
    self.findingIndexSliderWidget = ctk.ctkSliderWidget()
    self.findingIndexSliderWidget.singleStep = 1
    self.findingIndexSliderWidget.minimum = 0
    self.findingIndexSliderWidget.maximum = 10
    self.findingIndexSliderWidget.setDecimals(0)
    self.findingIndexSliderWidget.enabled = False
    self.findingIndexSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")

    parametersFormLayout.addRow("Finding Index: ", self.findingIndexSliderWidget)
    # zone line edit
    self.uidLineEdit = qt.QLineEdit()
    self.uidLineEdit.setReadOnly(True)
    parametersFormLayout.addRow("Zone:", self.uidLineEdit)

    # zone line edit
    self.zoneLineEdit = qt.QLineEdit()
    self.zoneLineEdit.setReadOnly(True)
    parametersFormLayout.addRow("Zone:", self.zoneLineEdit)


    # zone line edit
    self.clinSigLineEdit = qt.QLineEdit()
    self.clinSigLineEdit.setReadOnly(True)
    parametersFormLayout.addRow("Clinical Significance:", self.clinSigLineEdit)

    # connections
    self.loadButton.connect('clicked(bool)', self.onLoadButton)
    self.jsonFindingsPath.connect("currentPathChanged(QString)", self.onSelect)
    self.findingIndexSliderWidget.connect('valueChanged(double)', self.findingIndexChanged)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Load button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    path = self.jsonFindingsPath.currentPath
    self.loadButton.enabled = 'json' in path and os.path.isfile(path)

  def onLoadButton(self):
    import json
    with open(self.jsonFindingsPath.currentPath) as json_data:
            self.findings_json= json.load(json_data)
    self.findingIndexSliderWidget.maximum= len(self.findings_json)
    self.findingIndexSliderWidget.enabled = True
    self.findingIndexChanged(0)

  def load_file(self, file_path, node_name):
      result = slicer.util.loadVolume(file_path,returnNode=True)
      if result[0]:
          node = result[1]
          node.SetName(node_name)
      return node

  def findingIndexChanged(self, val):
    import CompareVolumes
    cvlogic = CompareVolumes.CompareVolumesLogic()

    scene = slicer.mrmlScene
    scene.Clear(0)
    lm = slicer.app.layoutManager()
    lm.setLayout(21)
    keys = self.findings_json.keys()
    finding_json = self.findings_json[keys[int(val)]]

    tra = finding_json['files']['t2_tse_tra']['filename'][0]
    sag = finding_json['files']['t2_tse_sag']['filename'][0]
    cor = finding_json['files']['t2_tse_cor']['filename'][0]
    adc = finding_json['files']['ADC']['filename'][0]
    bval = finding_json['files']['BVAL']['filename'][0]
    ktrans = finding_json['files']['Ktrans']['filename'][0]
    l, p, s = finding_json['lps']
    self.zoneLineEdit.setText(finding_json['zone'])
    self.uidLineEdit.setText(keys[int(val)])
    self.clinSigLineEdit.setText(finding_json['ClinSig'])
    markups = slicer.modules.markups
    markupsLogic = markups.logic()
    n = markupsLogic.AddFiducial(-l,-p,s)
    fidID = markupsLogic.GetActiveListID()
    fidNode = slicer.mrmlScene.GetNodeByID(fidID)
    # l.JumpSlicesToNthPointInMarkup(node.GetID(),0)

    adc_node = self.load_file(adc,'adc')
    tra_node = self.load_file(tra,'tra')
    sag_node = self.load_file(sag,'sag')
    cor_node = self.load_file(cor,'cor')
    adc_node = self.load_file(adc,'adc')
    bval_node = self.load_file(bval,'bval')
    ktrans_node = self.load_file(ktrans,'ktrans')

    redSliceWidget = lm.sliceWidget('Red')
    redSliceNode = redSliceWidget.mrmlSliceNode()
    compositeNode = redSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(tra_node.GetID())

    yellowSliceWidget = lm.sliceWidget('Yellow')
    yellowSliceNode = yellowSliceWidget.mrmlSliceNode()
    compositeNode = yellowSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(sag_node.GetID())

    greenSliceWidget = lm.sliceWidget('Green')
    greenSliceNode = greenSliceWidget.mrmlSliceNode()
    compositeNode = greenSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(cor_node.GetID())

    redpSliceWidget = lm.sliceWidget('Red+')
    redpSliceNode = redpSliceWidget.mrmlSliceNode()
    compositeNode = redpSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(adc_node.GetID())

    yellowpSliceWidget = lm.sliceWidget('Yellow+')
    yellowpSliceNode = yellowpSliceWidget.mrmlSliceNode()
    compositeNode = yellowpSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(bval_node.GetID())

    greenpSliceWidget = lm.sliceWidget('Green+')
    greenpSliceNode = greenpSliceWidget.mrmlSliceNode()
    compositeNode = greenpSliceWidget.mrmlSliceCompositeNode()
    compositeNode.SetBackgroundVolumeID(ktrans_node.GetID())

    redSliceNode.RotateToVolumePlane(tra_node)
    yellowSliceNode.RotateToVolumePlane(sag_node)
    greenSliceNode.RotateToVolumePlane(cor_node)
    redpSliceNode.RotateToVolumePlane(adc_node)
    yellowpSliceNode.RotateToVolumePlane(bval_node)
    greenpSliceNode.RotateToVolumePlane(ktrans_node)

    markupsLogic.JumpSlicesToNthPointInMarkup(fidNode.GetID(),n)

#
# FindingsScrollerLogic
#

class FindingsScrollerLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, findingIndex, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : findingIndex, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('FindingsScrollerTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


class FindingsScrollerTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_FindingsScroller1()

  def test_FindingsScroller1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = FindingsScrollerLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
