from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())
# measure process time
t0p = timer.clock()
# measure wall time
t0w = timer.time()

from math import sin
mypi = 3.14159265359	

#
def StartTimeMeasuring():
    # measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")
        
def ApplyLoad( model_part):
    a = 200.0
    b = 200.0
    for node in model_part.Nodes:      
        xpos = float(node.X)
        ypos = float(node.Y)
        value = -1e7*sin(mypi*xpos/a)*sin(mypi*ypos/b)
        node.SetSolutionStepValue(SURFACE_LOAD,0,[0,0,value])

#### TIME MONITORING END ####

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


######################################################################################
######################################################################################
######################################################################################

#### PARSING THE PARAMETERS ####

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#### model_part settings start ####

#defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

#add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
solver.AddVariables()

#read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

#add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# if we integrate it in the model part we cannot use combined solvers
solver.AddDofs()

#build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
##TODO: replace MODEL for the Kratos one ASAP
##get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#### model_part settings end ####


#### processes settings start ####

#obtain the list of the processes to be applied

import process_factory
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )

list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

#list_of_processes = []
#process_definition = ProjectParameters["boundary_conditions_process_list"]
#for i in range(process_definition.size()):
#    item = process_definition[i]
#    module = __import__(item["implemented_in_module"].GetString())
#    interface_file = __import__(item["implemented_in_file"].GetString())
#    p = interface_file.Factory(item, Model)
#    list_of_processes.append( p )
#    print("done ",i)
            
#print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

#TODO: decide which is the correct place to initialize the processes 
for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####

#### START SOLUTION ####

#TODO: think if there is a better way to do this
computing_model_part = solver.GetComputingModelPart()


#### output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()

# restart write included in gid IO ??

#### output settings end ####

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()

print(" ")
print("::[KSM Simulation]:: Analysis -START- ")

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
## Set results when are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

## Stepping and time settings (get from process info or solving info)
#delta time
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
#start step
step       = 0
#start time
time       = ProjectParameters["problem_data"]["start_time"].GetDouble()
#end time
end_time   = ProjectParameters["problem_data"]["end_time"].GetDouble()

# monitoring info:  # must be contained in the solver
#import solving_info_utility as solving_info_utils
#solving_info = solving_info_utils.SolvingInfoUtility(model_part)

# writing a initial state results file (if no restart)
# gid_io.write_results(time, computing_model_part) done in ExecuteBeforeSolutionLoop()
curvaturefile = open("curvatures.txt", "w")
momentfile = open("moments.txt", "w")
dispfile = open("displacements.txt", "w")
rotationfile = open("rotations.txt", "w")
stress_top_surface_file = open("stress_top_surface.txt", "w")
stress_bottom_surface_file = open("stress_bottom_surface.txt", "w")
# solving the problem (time integration)
while(time <= end_time):

    #TODO: this must be done by a solving_info utility in the solver
    # store previous time step
    #~ computing_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = delta_time
    # set new time step ( it can change when solve is called )
    #~ delta_time = computing_model_part.ProcessInfo[DELTA_TIME]

    time = time + delta_time
    step = step + 1
    main_model_part.CloneTimeStep(time)

    # print process info
    ##
    print("==================================================\nAPPLYING CUSTOM SIN SURFACE LOAD\n==================================================\n\n")
    ApplyLoad(main_model_part)
    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
        
    solver.Solve()
    
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    gid_output.ExecuteFinalizeSolutionStep()
    
        #TODO: decide if it shall be done only when output is processed or not (boundary_conditions_processes ??)
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # write results and restart files: (frequency writing is controlled internally by the gid_io)
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()
    
    # Displacements and rotations
    for node in main_model_part.Nodes:
        if node.Y > 99.9 and node.Y < 100.1:
            dispfile.write(str(node.X) + '\t' + str(node.GetSolutionStepValue(DISPLACEMENT_Z,0)) + "\n")
            rotationfile.write(str(node.X) + '\t' + str(node.GetSolutionStepValue(ROTATION_Y,0)) + "\n")

            
    # Gauss results
    proc_info = main_model_part.ProcessInfo 
    for element in main_model_part.Elements:
        x = 0.0
        y = 0.0
        for node in element.GetNodes():
            x += node.X/4.0
            y += (node.Y/4.0)
        if y > 90.0 and y < 100:    #only take elements along the top boundary
            # Curvatures
            curvaturefile.write(str(x) + "\t" + str(y))
            strain_result = []
            strain_result = element.GetValuesOnIntegrationPoints(SHELL_CURVATURE_GLOBAL, proc_info)
            strain_av = [0,0,0,0,0,0,0,0,0]
            for gauss_point in range(4):    #average gauss point values into central one
                for i in range (9):
                    strain_av[i] += strain_result[gauss_point][i]/4.0
            for i in range(9):
                curvaturefile.write("\t" + str(strain_av[i]))
            curvaturefile.write(("\n"))
            
            # Moments
            momentfile.write(str(x) + "\t" + str(y))
            strain_result = []
            strain_result = element.GetValuesOnIntegrationPoints(SHELL_MOMENT_GLOBAL, proc_info)
            strain_av = [0,0,0,0,0,0,0,0,0]
            for gauss_point in range(4):    #average gauss point values into central one
                for i in range (9):
                    strain_av[i] += strain_result[gauss_point][i]/4.0
            for i in range(9):
                momentfile.write("\t" + str(strain_av[i]))
            momentfile.write(("\n"))
            
            # Top surface stress
            stress_top_surface_file.write(str(x) + "\t" + str(y))
            strain_result = []
            strain_result = element.GetValuesOnIntegrationPoints(SHELL_STRESS_TOP_SURFACE_GLOBAL, proc_info)
            strain_av = [0,0,0,0,0,0,0,0,0]
            for gauss_point in range(4):    #average gauss point values into central one
                for i in range (9):
                    strain_av[i] += strain_result[gauss_point][i]/4.0
            for i in range(9):
                stress_top_surface_file.write("\t" + str(strain_av[i]))
            stress_top_surface_file.write(("\n"))
            
            # Bottom surface stress
            stress_bottom_surface_file.write(str(x) + "\t" + str(y))
            strain_result = []
            strain_result = element.GetValuesOnIntegrationPoints(SHELL_STRESS_BOTTOM_SURFACE_GLOBAL, proc_info)
            strain_av = [0,0,0,0,0,0,0,0,0]
            for gauss_point in range(4):    #average gauss point values into central one
                for i in range (9):
                    strain_av[i] += strain_result[gauss_point][i]/4.0
            for i in range(9):
                stress_bottom_surface_file.write("\t" + str(strain_av[i]))
            stress_bottom_surface_file.write(("\n"))
                      
    #TODO: decide if it shall be done only when output is processed or not
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


for process in list_of_processes:
    process.ExecuteFinalize()

curvaturefile.close()
dispfile.close()
momentfile.close()
rotationfile.close()
stress_top_surface_file.close()
# ending the problem (time integration finished)
gid_output.ExecuteFinalize()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

# check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.

#### END SOLUTION ####

# measure process time
tfp = timer.clock()
# measure wall time
tfw = timer.time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)






