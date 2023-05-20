import vtk
import os

# Define file paths
current = os.path.dirname(os.path.realpath(__file__))
vtk_loc = os.path.join(current, 'VTK_resampled', 'VTK_resampled', 'boone-pickens_')
out_loc = os.path.join(current, 'VTK_resampled')
file_extension = ".vtk"

# Read the input file
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName(input_file)
reader.Update()

# Create the output file
writer = vtk.vtkXdmfWriter()
writer.SetFileName(output_file)
writer.SetInputConnection(reader.GetOutputPort())
writer.Write()