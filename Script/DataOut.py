import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
np.set_printoptions(threshold=sys.maxsize)

class DataOut(object):
        def __init__(self): #No need to implement
             pass
         
        def _out2xyz_(self,label,position,filenamePre,tag):
            atomNum=np.shape(position)[1]
            
            for i in range(0,np.shape(position)[0]):
                filename = 'outxyz/'+filenamePre+tag+'/'+str(i)+'.xyz'
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                all_info=np.column_stack((label,position[i]))                
                print(str(atomNum)+'\n',file = open(filename,"a"))

                for k in range(0,atomNum):
                    print (*all_info[k], file = open(filename,"a")) 
                    
            return
        
        
        def _out2cif_(self,fileFolder, label, cellParameter, position, charge, filenamePre):
            # mkdir
            if charge=='none':
                charge=np.zeros([np.shape(label)[0],1])                
                label=np.array(label).reshape([np.shape(label)[0],1])
            if np.shape(position)[1]==3:
                positionT=position[:].copy()
                position=np.expand_dims(positionT, axis=0)
            
            for i in range(0,np.shape(position)[0]):
                #if np.shape(position)[2]==3:
                    #filename = fileFolder+'/'+'RigidOrSingleOut'+'.cif'
                #else:
                filename = fileFolder+'/'+filenamePre+'/'+str(i)+'.cif'
                    
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                all_info=np.column_stack((label,label,position[i],charge))
            #           
                t=list(['_cell_length_a '+str(cellParameter[0])])
                t.append('_cell_length_b '+str(cellParameter[1]))
                t.append('_cell_length_c '+str(cellParameter[2]))
                t.append('_cell_angle_alpha '+str(cellParameter[3]))
                t.append('_cell_angle_beta '+str(cellParameter[4]))
                t.append('_cell_angle_gamma '+str(cellParameter[5]))
                t.append("_symmetry_space_group_name_H-M         'P 1'")
                t.append('_symmetry_Int_Tables_number            1')
                t.append('loop_')
                t.append('_symmetry_equiv_pos_as_xyz')
                t.append("   'x, y, z'")
                t.append('                         ')
                t.append('loop_')
                t.append('_atom_site_label') 
                t.append('_atom_site_type_symbol') 
                t.append('_atom_site_fract_x ')
                t.append('_atom_site_fract_y')
                t.append('_atom_site_fract_z')
                #t.append('_atom_site_occupancy')
                t.append('_atom_site_charge')
                        
                position=np.array(position)
                for j in range (0,19):
                    print(t[j],file = open(filename,"a"))
                    
                for k in range(0,np.shape(position)[1]):
                    print (*all_info[k,:], file = open(filename,"a")) 
        
            return
        
        def _outParityPlot4cmp_(snapAtomName,cifBasedVib,chargeBasedstd,mofName):
            chargeBasedstd_iso=[]
            nameAxis=[]
            cifAxis=[]
            for i in range(0,len(snapAtomName)):
                for j in range(0,len(chargeBasedstd[i])):
                    temp=np.sqrt(np.mean(chargeBasedstd[i][j,:]**2))
                    nameAxis.append(snapAtomName[i].copy())
                    chargeBasedstd_iso.append(temp)
            tagSnap=np.chararray((len(nameAxis), 1))
            tagSnap[:]='Snap'
            dataSnap=np.column_stack((np.array(nameAxis,dtype=str),np.array(chargeBasedstd_iso,dtype=float),tagSnap))
            
            for i in range(0,len(cifBasedVib)):
                cifAxis.append(cifBasedVib[i,0][0])
            cifBasedVib_iso=(np.sqrt(np.array(cifBasedVib[:,1],dtype=float)))
            tagCif=np.chararray((len(cifAxis), 1))
            tagCif[:]='cif'
            dataCif=np.column_stack((np.array(cifAxis,dtype=str),np.array(cifBasedVib_iso,dtype=float),tagCif))
            dataAll=np.concatenate((dataSnap,dataCif))
            
            dataIn=pd.DataFrame(dataAll,columns=['name','Vib_iso','tag'])
            dataIn['Vib_iso']=dataIn['Vib_iso'].apply(pd.to_numeric)
            #fig, ax= plt.subplots(figsize=(6, 6))
            sns.set(style="ticks")
            sns.catplot(x='name', y='Vib_iso',hue="tag",data=dataIn)
            
            plt.savefig('figure/cmp'+mofName+'.jpg',dpi=500)
            
            return 
            
            
            