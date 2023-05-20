import h5py as h5
import numpy as np
import configparser
import os
import json
import matplotlib.pyplot as plt
import sys
import pandas as pd
import openpyxl

from pathlib import Path
# sys.path is a list of absolute path strings
# sys.path.append(str(Path.home()) + '\Desktop\python\sharpy')

from OSU_Contribution.FEM import FEM
from OSU_Contribution.Aero import Aero

abss = os.path.dirname(os.path.abspath(__file__))

class Paramaterized():
    def __init__(self, case_name, route, file, xfoil, excel = False):
        self.case_name = case_name
        self.file = file
        self.xfoil = xfoil
        self.excel = excel
        self.route = route

        self.clean_test_files(self.case_name, self.route)
        
        if self.excel == True:
            data = self.file 
        
            # Creates FEM HF5 File to pass the variables that are needed in the aero creator 
            pass_aero = FEM(data, route, self.case_name, self.excel)
            self.app_nodes = pass_aero.app_nodes # to be used for static trim settings 
            self.ref_lookup = pass_aero.ref_lookup
            self.lumped_loc = pass_aero.lumped_loc
            self.beam_opt_nodes = pass_aero.beam_opt_nodes
            self.beam_opt_elems = pass_aero.beam_opt_elems
            self.beam_opt_sys = pass_aero.beam_opt_sym
            self.num_nodes = pass_aero.num_node
            
            df = pd.read_excel(data, sheet_name= 'Airfoils')
            rows = df.shape[0]
            columns = df.shape[1]
            
            air_data = []
            for i in range(0, columns, 2):
                air_temp = []
                for j in range(1, rows):
                    if str(df.values[j, i]) != 'nan':
                        air_temp.append([df.values[j, i], df.values[j, i + 1]])
                    if str(df.values[j, i]) == 'nan' and air_temp != []:
                        air_data.append(air_temp)
                        break
                    if j == rows - 1 and air_temp != []:
                        air_data.append(air_temp)
            
            pass_run = Aero(route, self.case_name, data, air_data, pass_aero, xfoil, self.excel)
            self.tail_index = pass_run.tail_index 
        
        if self.excel == False:
            f = open(self.file)
            data = json.load(f)
            
            # Creates FEM HF5 File to pass the variables that are needed in the aero creator 
            pass_aero = FEM(data, route, self.case_name, self.excel)
            
            # loading airfoils from 
            ah95 = self.airfoil_reader(abss + '\\ah95160c.txt')
            MH212 = self.airfoil_reader(abss + '\\MH212_2.txt')
            NACA63 = self.airfoil_reader(abss + '\\NACA63015F.txt')

            air_data = []

            air_data.append(ah95)
            air_data.append(MH212)
            air_data.append(NACA63)
            
            # Creates Aero HF5 File
            Aero(route, self.case_name, data, air_data, pass_aero, xfoil, self.excel)


    def airfoil_reader(self, Name):
        airfoil = []
        f = open(Name, "r")
        for line in f:
            current = line.strip('\n')
            current = current.split(" ")
            final = []
            for i in range(len(current)):
                if current[i] != '':
                    final.append(current[i])
            current = final
            try:
                one = float(current[0])
                two = float(current[1])
                if type(float(current[1])) == float:
                    airfoil.append([one, two])
            except:
                pass
        return np.array(airfoil)

    def clean_test_files(self, case_name, route):
        fem_file_name = os.path.join(route, '{}.fem.h5'.format(case_name))
        if os.path.isfile(fem_file_name):
            os.remove(fem_file_name)

        solver_file_name = os.path.join(route, '{}.sharpy'.format(case_name))
        if os.path.isfile(solver_file_name):
            os.remove(solver_file_name)
        
        aero_file_name = os.path.join(route, '{}.aero.h5'.format(case_name))
        if os.path.isfile(aero_file_name):
            os.remove(aero_file_name)

        dyn_file_name = os.path.join(route, '{}.dyn.h5'.format(case_name))
        if os.path.isfile(dyn_file_name):
            os.remove(dyn_file_name)
            
            
    def excel_write(excel_og, new_excel, paramater, section, new_value):
        """Takes input file template and make a new input file with updated values

        Args:
            excel_og (str): Name of the excel file template. MAKE SURE that the file is in same directory as Run_Sarpy.py
                                                Ex: 'test.xlsx'. DO NOT INCLUDE PATH
            new_excel (str): Name of the new excel file. Ex: 'New_File.xlsx'
            paramater (list[str]): The name of the paramater that will be changed. Parameter name is the same as
                                                the title of the column in excel file template
            section (list[int]): the excel rows which values will be changed. Get the row numbers from the 
                                                excel sheet template
            new_value (list[float]): new values for the paramater of choice
        """
        direct = __file__.split('/')
        del direct[-1]
        direct = '/'.join(direct)
        # new file creation 
        # validates that there is no directory then creates one
        path = os.path.join(direct + '/',  "Created_Input_Files")
        if os.path.isdir(path) != True:
            os.mkdir(path)
            
        # checks to see if file of the same name is already there
        if os.path.isfile(direct + '/Created_Input_Files' + '/' + excel_og):
            os.remove(direct + '/Created_Input_Files' + '/' + excel_og)
        
        # copies reference excel into new directory
        import shutil  
        shutil.copy(direct + '/' + excel_og, direct + '/Created_Input_Files')
        
        # checks to see if file of the same name is already there and removes
        if os.path.isfile(direct + '/Created_Input_Files' + '/' + new_excel):
            os.remove(direct + '/Created_Input_Files' + '/' + new_excel)
        
        # renaming copied excel template
        os.rename(direct + '/Created_Input_Files' + '/' + excel_og, direct + '/Created_Input_Files' + '/' + new_excel)

        excel_og = direct + '/' + excel_og
        master = pd.read_excel(excel_og, sheet_name= 'Master')
        master_titles = master.values[0, 1:]
        #removing units parenthisis from title
        for i in range(len(master_titles)):
            temp = str(master_titles[i]).lower().split(' ')
            master_titles[i] = temp[0]
        
        # defining dictionary 
        dict_master = {} 
        dict_master_index = {}
        # capturing all of the values in master
        for i in range(master.shape[0]): 
            for j in range(1, master.shape[1]-1):
                if i == 0:
                    dict_master[master_titles[j]] = []
                    dict_master_index[master_titles[j]] = j + 2
                else:
                    dict_master[master_titles[j]].append(master.values[i, j])
                    
                    
        lumped = pd.read_excel(excel_og, sheet_name= 'Applied_Forces_Lumped_Masses')
        lumped_titles = lumped.values[0, 0:9]
        lumped_titles = np.append(lumped_titles, lumped.values[0, 11:])
        
        dict_lumped = {}
        dict_lumped_index = {}
        #removing units parenthisis from title
        for i in range(len(lumped_titles)):
            temp = str(lumped_titles[i]).lower().split(' ')
            lumped_titles[i] = str(temp[0]).lower()
            
        # capturing all of the values in lumped
        for i in range(lumped.shape[0]): 
            for j in range(0, 9):
                if i == 0:
                    dict_lumped[lumped_titles[j]] = []
                    dict_lumped_index[lumped_titles[j]] = j
                else:
                    # does not append nan to dictionary
                    if str(lumped.values[i, j]) == 'nan':
                        continue
                    else:
                        dict_lumped[lumped_titles[j]].append(lumped.values[i, j])
                    
            for j in range(10, lumped.shape[1]):
                if i == 0:
                    dict_lumped[lumped_titles[j-2]] = []
                    dict_lumped_index[lumped_titles[j-2]] = j
                else:
                    # does not append nan to dictionary
                    if str(lumped.values[i, j]) == 'nan':
                        continue
                    else:
                        dict_lumped[lumped_titles[j-2]].append(lumped.values[i, j])

        control = pd.read_excel(excel_og, sheet_name= 'Control_Surfaces')
        control_titles = control.values[0, 0:]
        
        #removing units parenthisis from title
        for i in range(len(control_titles)):
            temp = str(control_titles[i]).lower().split(' ')
            control_titles[i] = temp[0]
        
        dict_control = {}
        dict_control_index = {}
        # capturing all of the values in control
        for i in range(control.shape[0]): 
            for j in range(control.shape[1]):
                if i == 0:
                    dict_control[control_titles[j]] = []
                    dict_control_index[control_titles[j]] = j
                else:
                    dict_control[control_titles[j]].append(control.values[i, j])
                    
        dict_update = {}
        dict_update_index = {}
        update = pd.read_excel(excel_og, sheet_name= 'Update')
        update_titles = update.values[0, 1:]

        #removing units parenthisis from title
        for i in range(len(update_titles)):
            temp = str(update_titles[i]).lower().split(' ')
            update_titles[i] = temp[0]
        
        stop = 0
        for i in range(update.shape[0]): 
            if stop == 1:
                break
            for j in range(0, update.shape[1]-1):
                if i == 0:
                    dict_update[update_titles[j]] = []
                    dict_update_index[update_titles[j]] = j + 2
                else:
                    if str(update.values[i,1]).lower() == 'nan':
                        stop = 1
                        break
                    else:
                        dict_update[update_titles[j]].append(update.values[i, j+1])
            

        excel_wb = openpyxl.load_workbook(direct + '/Created_Input_Files' + '/' + new_excel)
        test = excel_wb.worksheets[0]
        # writing new values into dictionary
        for i in paramater:
            i = str(i).lower()
            for j in range(len(section)):
                for k in range(len(section[0])):
                # -3 due to being in reference to excel sheet
                    # finds the correct sheet to update to 
                    try:
                        dict_control[i][int(section[j][k]-3)] = new_value[j][k]
                        sheet = 3
                        dic = dict_control
                        dic_index = dict_control_index
            
                    except KeyError:    
                        try:
                            dict_lumped[i][int(section[j][k]-3)] = new_value[j][k]
                            sheet = 2
                            dic = dict_lumped
                            dic_index = dict_lumped_index
            
                        except KeyError:
                            try:
                                dict_update[i][int(section[j][k]-3)] = new_value[j][k]
                                sheet = 0
                                dic = dict_update
                                dic_index = dict_update_index
                                
                            except KeyError:
                                sheet = 5
                                dic = dict_master
                                dic_index = dict_master_index
                                
                    # updating the excel sheet
                    excel_wb.worksheets[sheet].cell(row = int(section[j][k]), column = dic_index[i]).value = dic[i][k]
                    
        # saves the excel sheet with the updated values  
        excel_wb.save(direct + '/Created_Input_Files' + '/' + new_excel)  


