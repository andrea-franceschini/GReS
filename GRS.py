import tkinter as tk
import os 
from tkmacosx import Button as bttn
import gmsh

import faulthandler
faulthandler.enable()

path_to_simParam = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/simParam.dat'
path_to_materialsList2 = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/materialsList2.dat'
path_to_dirBCSurf = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/dir_BCSurf_poro_tetra.dat'
path_to_mod_dirBCSurf = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/modified_dir_BCSurf_poro_tetra'
path_to_neuSurf = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/neuSurf_BC_poro_tetra.dat'
path_to_mod_neuSurf = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/modified_neuSurf_BC_poro_tetra'
path_to_main = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp/main.m'

gmsh.initialize()

def raises(name):
    name.tkraise()
    match str(name):
        case ".!frame2.!frame":
            config_simParam()
        case ".!frame2.!frame2":
            config_Materials()
        case ".!frame2.!frame3":
            config_conditions()
        case ".!frame2.!frame4":
            config_Models()
            
def frame_destroyer():
    for widget in simParam_frame.winfo_children():
        widget.destroy()
    for widget in times_frame.winfo_children():
        widget.destroy()
    for widget in models_frame.winfo_children():
        widget.destroy()
    for widget in materials_frame.winfo_children():
        widget.destroy()          

####################################### SIM PARAMETERS ###################################################
def config_simParam():
    #for widget in simParam_frame.winfo_children():
    #    widget.destroy()
    frame_destroyer()
    print("opened")
    apply_changes_frame = tk.Frame(simParam_frame)
    apply_changes_frame.grid(row=15, column=0, columnspan=3)
    #Per modificare sim_param: basta leggere il contenuto del file e modificarlo
    #filename = 'simParam.dat'
    #path = '/Users/alessandrodoro/Downloads/GReS/GReS/Tests/mechanics_dp'
    file = open(path_to_simParam, 'r+')
    i = 0
    text_list = []
    label_list = []
    desc_list = []
    for line in file:
        if line != 'End':
            value, desc = line.split(' ', 1)
            print(value, desc)
            desc_list.append(desc)
            text_list.append(tk.Entry(simParam_frame))
            text_list[i].insert(0, value)
            text_list[i].grid(row=i, column=0)
            label_list.append(tk.Label(simParam_frame, text=desc))
            label_list[i].grid(row=i, column=1)
            i += 1
        else:
            pass
  
    def rewrite(i):
        j = 0
        file = open(path_to_simParam, 'w+')
        for j in range(i):
            newvalue = text_list[j].get()
            newline = str(newvalue)+' '+str(desc_list[j])
            file.write(newline)
        file.write("End")
        file.close()
        
    apply_button = bttn(apply_changes_frame, text='Apply changes', command=lambda: rewrite(i), bg='#567EE7', fg='white', activebackground='cyan', borderless=1)
    apply_button.grid(row=1, column=1, columnspan=3)
    
    
    
############################################ MATERIALS ######################################################
def config_Materials():
    frame_destroyer()
    file = open(path_to_materialsList2, 'r+')
    name = 'PorousMedia2.dat' #Elastic Class
    name1 = 'PorousMedia3.dat' #Druker Prager Class
    for line in file:
        if line != 'End':
          desc = line.split('/')
        else:
           pass
    desc[1] = desc[1].replace("\n","")
    print(desc)
    if desc[1] == name:
        label = tk.Label(materials_frame, text="You're in Elastic Class")
    else:
        label = tk.Label(materials_frame, text="You're in Drucker Prager Class")
    label.grid(row=0, column=1, columnspan=3)
    
    def check():
        file = open(path_to_materialsList2, 'r+')
        if not file.read():
            print('file vuoto')
        file.seek(0)
        for line in file:
            if line != 'End':
                desc = line.split('/')
            else:
                pass
        desc[1] = desc[1].replace("\n","")
        if desc[1] == name:
            newline = '1 Materials/'+ name1 + '\n'
            label.config(text="You're in Drucker Prager Class")
        else:
            newline = '1 Materials/'+ name + '\n'
            label.config(text="You're in Elastic Class")
        
        file = open(path_to_materialsList2, 'w')
        file.write(newline)
        
        file.write('End')
        file.close()
    
    check_button = bttn(materials_frame, text="Change class", command = check, borderless=1)
    check_button.grid(row=3, column=1, columnspan=3)
################################## CONDITIONS ############################################
def config_conditions():
    frame_destroyer()
    for widget in times_frame.winfo_children():
        widget.destroy()
    dir_button = tk.Button(times_frame, text='Dirichlet Conditions', command=lambda: config_Times('dir'))
    dir_button.grid(row=0, column=0)
    neu_button = tk.Button(times_frame, text='Neumann Conditions', command=lambda: config_Times('neu'))
    neu_button.grid(row=0, column=1)
    
    
    

######################################## TIMES #####################################################
def config_Times(path_to):
    frame_destroyer()
    #for widget in times_frame.winfo_children():
    #    widget.destroy()
    
    if path_to == 'dir':
        path_to_dir = path_to_dirBCSurf
        path_to_modified_dir = path_to_mod_dirBCSurf
        mod_path = 'modified_dir_BCSurf_poro_tetra'
        k = 5
    else: 
        path_to_dir = path_to_neuSurf
        path_to_modified_dir = path_to_mod_neuSurf
        mod_path = 'modified_neuSurf_BC_poro_tetra'
        k = 6

    main_dir_file = open(path_to_dir, 'r')
    
    def remove(button_id):
        os.remove(path_to_modified_dir+'/time'+str(button_id)+'.dat')
        main_dir_file = open(path_to_dir, 'r')
        main_dir_file.seek(0)
        remove_line = str(float(button_id))+' '+mod_path+'/time'+str(button_id)+'.dat'
        print(remove_line)
        lines = main_dir_file.readlines()
        filtered_lines = [line for line in lines if remove_line not in line]
        print(filtered_lines)
        main_dir_file = open(path_to_dir, 'w')
        main_dir_file.writelines(filtered_lines)
        print('file deleted')
        for item in times_frame.winfo_children():
            if item.grid_info().get('row') == button_id:
                item.destroy()
        
    
    def modify(button_id):
        def get_value(newvalue):
            list_file = open(path_to_modified_dir+'/list','r')
            list_file.seek(0)
            list_size = len(list_file.readlines())
            time_file = open(path_to_modified_dir+'/time'+str(button_id)+'.dat','r')
            time_file.seek(0)
            all_time_lines = time_file.readlines()
            first_time_line = all_time_lines[0]
            newlines = []
            for i in range(list_size-2):
                newlines.append(str(newvalue.get())+'\n')
            newlines.append(str(newvalue.get()))
            newlines.insert(0, first_time_line)
            time_file = open(path_to_modified_dir+'/time'+str(button_id)+'.dat','w')
            time_file.writelines(newlines)
            time_file.close()
            list_file.close()
            for widget in times_frame.winfo_children():
                if widget.grid_info().get('row') == button_id:
                    if widget.grid_info().get('column') == 4 or widget.grid_info().get('column') == 5:
                        widget.destroy()

                    
        newvalue = tk.Entry(times_frame, text=inside_stuff[1])
        newvalue.grid(row=button_id, column=4)
        get_newvalue = bttn(times_frame, text="Apply", command = lambda: get_value(newvalue), fg = 'White', bg = '#567EE7')
        get_newvalue.grid(row=button_id, column=5)
        print(button_id)
    
    def add(i):
        
        ### Add file ###
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        params = open(path_to_simParam, 'r')
        params.seek(0)
        params_lines = params.readlines()
        list_params_lines = []
        for lines in params_lines:
            lines = lines.split(' ',1)
            list_params_lines.append(lines)
        end_time = int(float(list_params_lines[0][0]))
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        add_button.destroy()
        new_value = tk.Scale(times_frame, activebackground='#cdf2f7', orient='horizontal', from_=-100, to=100, sliderlength=30, width=10, length=180)
        new_value.grid(row=i+1, column=1, columnspan=3)
        newvalue_label = tk.Label(times_frame, text='Force: ')
        newvalue_label.grid(row=i+1, column=0)
        new_time = tk.Scale(times_frame, activebackground='#cdf2f7', orient='horizontal', from_=0, to=end_time, sliderlength=30, width=10, length=180)
        new_time.grid(row=i+2, column=1, columnspan=3)
        newtimelabel = tk.Label(times_frame, text='Time: ')
        newtimelabel.grid(row=i+2, column=0)
        list_file = open(path_to_modified_dir+'/list','r')
        list_file.seek(0)
        list_size = len(list_file.readlines())
        def terminal_add():
            time_file = open(path_to_modified_dir+'/time'+str(new_time.get())+'.dat','w')
            firstline = '% TIME '+str(new_time.get())+'.000000'+'\n'
            newlines = []
            for j in range(list_size-2):
                newlines.append(str(new_value.get())+'\n')
            newlines.append(str(new_value.get()))
            newlines.insert(0, firstline)
            time_file.writelines(newlines)
            time_file.close()
            ### Add Reference in file ### 
            main_dir_file = open(path_to_dir, 'r')
            time_lines = main_dir_file.readlines()
            first_time_lines = time_lines[0:k]
            time_lines = time_lines[k:-1]
            time_dict = {}
            for line in time_lines:
                line = line.split(' ', 1)
                time_dict[float(line[0])] = line[1]
            time_dict[float(new_time.get())] = mod_path+'/time'+str(new_time.get())+'.dat'+'\n'
            time_dict = dict(sorted(time_dict.items(), key=lambda x:x[0]))
            final_lines = []
            for key, value in time_dict.items():
                line = str(key)+' '+value
                final_lines.append(line)
            final_lines.append('End')
            final_lines = first_time_lines+final_lines
            main_dir_file = open(path_to_dir, 'w')
            main_dir_file.writelines(final_lines)
            main_dir_file.close()
            added_button.destroy()
            new_value.destroy()
            newvalue_label.destroy()
            new_time.destroy()
            newtimelabel.destroy()
            config_Times(path_to)
            
        added_button = bttn(times_frame, text='Add Time', bg='#567EE7', fg='white', command = terminal_add)
        added_button.grid(row=i+3, column=1, columnspan=5)
        
    content = main_dir_file.readlines()
    remove_buttons_list = []
    modify_buttons_list = []
    times_labels = []
    inside_stuff_labels = []
    i = 0
    
    for line in content[k:-1]:
        value = line.split('/')
        time_list = value[0].split(" ")
        time_list[0] = float(time_list[0])
        time = int(time_list[0])
        #print('time', time)
        value[1] = value[1].replace('\n','')
        
        times = open(path_to_modified_dir+'/'+value[1], 'r')
        inside_stuff = times.readlines()
        
        times_labels.append(tk.Label(times_frame, text = value[1]))
        times_labels[i].grid(row=time, column=0)
        
        inside_stuff_labels.append(tk.Label(times_frame, text =str(inside_stuff[1])))
        inside_stuff_labels[i].grid(row=time, column=1)
        
        remove_buttons_list.append(bttn(times_frame, text='Remove', fg='white', bg='red', command=lambda r=time: remove(r)))
        remove_buttons_list[i].grid(row=time, column=2)
        
        modify_buttons_list.append(bttn(times_frame, text='Modify', fg='White', bg='Green', command=lambda c=time: modify(c)))
        modify_buttons_list[i].grid(row=time, column=3)
        
        i += 1
    add_button = bttn(times_frame, text='Add Time', command =lambda: add(time), fg='White', bg='#567EE7')
    add_button.grid(row=time+1, column=1, columnspan=5)
    
        
    
def config_Models():
    frame_destroyer()
    main_dir = '/Users/alessandrodoro/Downloads/TERZAGHI_MESHES'
    def find_file(main_dir):
        file_msh = []
        name_files = []
        for root, _, files in os.walk(main_dir):  
            for file in files:
                if file.endswith('.msh'):
                    file_msh.append(os.path.join(root, file))
                    name_files.append(file)           
        return file_msh, name_files
    
    def openandreplace(button_id):
        main_main_path = path_to_main
        id = str(button_id).split('/')
        print(id[-1], type(id[-1]))
        main_main = open(main_main_path, 'r')
        main_main.seek(0)
        main_lines = main_main.readlines()
        print('main_lines', main_lines)
        newline = "fileName = "+f"'{id[-1]}'"+";\n"
        main_lines[16] = newline
        main_main = open(main_main_path, 'w')
        main_main.writelines(main_lines)
        main_main.close()
        
        gmsh.open(str(button_id))
        physical_groups = gmsh.model.getPhysicalGroups()
        elementi_2d_per_gruppo = {}
        for dim, tag in physical_groups:
            if dim != 2:
                continue          
            name = gmsh.model.getPhysicalName(dim, tag)  
            entity_tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            elementi_2d_per_gruppo[name] = []
            for entity_tag in entity_tags:
                element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)    
                for elem_type, elems in zip(element_types, element_tags):
                    if elem_type == 2 or elem_type == 3:  # 2=Triangoli, 3=Quadrilateri
                        elementi_2d_per_gruppo[name].extend(elems)
        zipped_values = []
        for gruppo, elementi in elementi_2d_per_gruppo.items():
            if gruppo == 'bottom' or gruppo =='lateral':
                zipped_values.extend(elementi)
        num1 = len(zipped_values)
        num2 = len(elementi_2d_per_gruppo['bottom'])
        zipped_values.extend(zipped_values)
        zipped_values.extend(elementi_2d_per_gruppo['bottom'])
        first_line = str(num1)+' '+str(num1)+' '+str(num2)+"\t%#ID_constrained"+"\n"
        mod_zipped_values = [str(num)+'\n' for num in zipped_values[:-1]]
        mod_zipped_values.append(str(zipped_values[-1]))
        mod_zipped_values.insert(0, first_line)
        # questo finisce in modified_dir ^^^
        path = path_to_mod_dirBCSurf+'/list'
        os.remove(path)
        newlist = open(path, 'w')
        newlist.writelines(mod_zipped_values)
        print(type(elementi_2d_per_gruppo['top']), elementi_2d_per_gruppo['top'])
        num3 = len(elementi_2d_per_gruppo['top'])
        mod_zipped_values = [str(num)+'\n' for num in elementi_2d_per_gruppo['top'][:-1]]
        mod_zipped_values.append(str(elementi_2d_per_gruppo['top'][-1]))
        first_line = str(num3)+"\t%#ID_constrained"+"\n"
        mod_zipped_values.insert(0, first_line)
        newlist.close()
        # questo finisce in modified_neu ^^^
        path = path_to_mod_neuSurf+'/list'
        os.remove(path)
        newlist = open(path, 'w')
        newlist.writelines(mod_zipped_values)
        newlist.close()
        #-------------------- UPDATE TIMES -------------------------
        paths = [path_to_mod_dirBCSurf, path_to_mod_neuSurf]
        for path in paths:
            for root, _, files in os.walk(path):
                for file in files:
                    if file.endswith('.dat'):
                        name = os.path.join(root, file)
                        file = open(name, 'r') 
                        file.seek(0)
                        read_lines = file.readlines()
                        firstline = read_lines[0]
                        value = read_lines[1]
                        path_to_list = path+'/list'
                        list_file = open(path_to_list, 'r')           
                        list_file.seek(0)
                        list_files = list_file.readlines()
                        lenght_list = len(list_files)
                        final_lines = []
                        for i in range(lenght_list-2):
                            final_lines.append(value)
                        final_lines.append(value.strip('\n'))
                        final_lines.insert(0, firstline)
                        file.close()
                        file = open(name, 'w')
                        file.writelines(final_lines)
                        file.close()
                        list_file.close()
                        
    found_files, name_files = find_file(main_dir)
    i = 0
    file_label = []
    change_to = []
    for file in found_files:
        file_label.append(tk.Label(models_frame, text=name_files[i]))
        file_label[i].grid(row=i, column=0)
        change_to.append(bttn(models_frame, text='Change model', command=lambda i=file: openandreplace(i), bg='#567EE7', fg='white'))
        change_to[i].grid(row=i, column=1)
        i += 1
        
window = tk.Tk()

modes_frame = tk.Frame(window)
modes_frame.grid(row = 0, column = 0)
script_frame = tk.Frame(window)
script_frame.grid(row = 1, column = 0)

simParam_frame = tk.Frame(script_frame)
simParam_frame.grid(row=0, column=0, sticky='n')

materials_frame = tk.Frame(script_frame)
materials_frame.grid(row=0, column=0, sticky='')

times_frame = tk.Frame(script_frame)
times_frame.grid(row=0, column=0, sticky='n')

models_frame = tk.Frame(script_frame)
models_frame.grid(row=0, column=0, sticky='n')

simParam_button = bttn(modes_frame, text='SimParam', command=lambda: raises(simParam_frame), bg='#434343', fg='white')
simParam_button.grid(row=0, column=0)

materials_button = bttn(modes_frame, text='Materials', command=lambda: raises(materials_frame), bg='#434343', fg='white')
materials_button.grid(row=0, column=1)

times_button = bttn(modes_frame, text='Times', command=lambda: raises(times_frame), bg='#434343', fg='white')
times_button.grid(row=0, column=2)

models_button = bttn(modes_frame, text='Models', command=lambda: raises(models_frame), bg='#434343', fg='white')
models_button.grid(row=0, column=3)

window.geometry('610x600')
window.title('GReS Modifier')
window.mainloop()