from Script.CIFIO import CIF
import pymatgen as mg
import numpy as np
import re
import pandas as pd
from scipy.optimize import minimize
from Script.DataOut import DataOut
from scipy.spatial.transform import Rotation as R
from scipy.spatial import distance
import scipy

class cifOPandIO(object):
        def __init__(self,cifName,T_set): #No need to implement
            self.cifName = cifName
            self.T_set=T_set
            if cifName!=None:
                self._read_()
                print('start_reading_cif')
            
        def _read_(self,LigandRigid=True):
            """
            read clean
            """
            self.structure_clean = mg.Structure.from_file('exp_cifs/'+self.cifName+'_clean.cif')
            self.position_clean=self.structure_clean.cart_coords[:]
            a=self.structure_clean.lattice
            """
            read exp
            """            
            cifobj = CIF()
            cifobj.read('exp_cifs/'+self.cifName+'.cif')
            self.data = cifobj._data
            self.T_ratio=np.sqrt(self.T_set/float(self.data['_cell_measurement_temperature']))
            self.positionFromCif,self.labelFromCif,self.atomTypeFromCif,self.patternFromCif =self._symOp_()
            self.index_associatedH,self.index_H=self._findHList_(self.positionFromCif,self.atomTypeFromCif)
            
        def stdOutStruct(self):    
            self.gaussianStruct,b=self._fullRank_add_Gaussian_(self._frac2Cartesian_(self.positionFromCif),self.patternFromCif,self.index_associatedH,self.index_H,self.transM,10)
            self._out2cif_(self.gaussianStruct,self.atomTypeFromCif,'Gaussian')
            return
        
        def ellipsoidOut(self):     
            self.pattern_aniso=np.array(np.column_stack([self.data['_atom_site_aniso_U_11'],self.data['_atom_site_aniso_U_22'],self.data['_atom_site_aniso_U_33'],self.data['_atom_site_aniso_U_12'],self.data['_atom_site_aniso_U_13'],self.data['_atom_site_aniso_U_23']]),dtype=float)
            self.pattern_aniso_label=np.array([x.strip() for x in self.data['_atom_site_aniso_label']],dtype=str)
            self.gaussianStruct_elipsoid,b=self._fullRank_add_Gaussian_aniso(self._frac2Cartesian_(self.positionFromCif),self.index_associatedH,self.index_H,self.transM,10)
            self._out2cif_(self.gaussianStruct_elipsoid,self.atomTypeFromCif,'ellipsoid')
            
            return
        
        def CoREOutStruct(self):   
            cifobj = CIF()
            cifobj.read('exp_cifs/'+self.cifName+'_clean2.cif')
            label=np.array([re.split(r'[;,\s]\s*', x)[0] for i,x in enumerate(cifobj._data['_atom_site_label'])],dtype=str)
            CoREPattern=np.zeros((np.shape(label)[0],3))
            for i in range(0,np.shape(label)[0]):
                mark=np.where(label[i]==self.labelFromCif)
                if len(mark[0])!=0:
                    CoREPattern[i,:]=self.patternFromCif[mark[0][0]]
                
            self.index_associatedH_CoRE,self.index_H_CoRE=self._findHList_(self.position_clean,np.array(self.structure_clean.species,dtype=str))
            self.CoREgaussianStruct,b=self._fullRank_add_Gaussian_(self.position_clean,CoREPattern,self.index_associatedH_CoRE,self.index_H_CoRE,self.transMCoRE,50)
            self._out2cif_(self.CoREgaussianStruct,np.array(self.structure_clean.species,dtype=str),'CoRE',self.cellParameterCoRE)
            return 
        
        def matchStruct(self):
            if np.shape(self.positionFromCif)[0] >= np.shape(self.position_clean)[0]:
                labled=self._matchStructure_(self.positionFromCif,self.labelFromCif,self.atomTypeFromCif,self.structure_clean)
            else:
                raise "exp_cif is too small"
                
            return labled
            

        def _symOp_(self):
            sym=pd.DataFrame([re.split(r'[;,\s]\s*', x)[0:3] for i,x in enumerate(self.data['_symmetry_equiv_pos_as_xyz'])],columns=['x', 'y', 'z'])
            x=np.array(self.data['_atom_site_fract_x'],dtype=float)
            y=np.array(self.data['_atom_site_fract_y'],dtype=float)
            z=np.array(self.data['_atom_site_fract_z'],dtype=float)
            position=np.column_stack([x,y,z])
            label=np.array([re.split(r'[;,\s]\s*', x)[0] for i,x in enumerate(self.data['_atom_site_label'])],dtype=str)
            atomType=np.array([re.split(r'[;,\s]\s*', x)[0] for i,x in enumerate(self.data['_atom_site_type_symbol'])],dtype=str)
            pattern=np.sqrt(np.array(np.column_stack([self.data['_atom_site_U_iso_or_equiv'],self.data['_atom_site_U_iso_or_equiv'],self.data['_atom_site_U_iso_or_equiv']]),dtype=float))
            index_H=np.where( atomType == 'H')
            labelOrigin=label.copy()
            atomTypeOrigin=atomType.copy()
            patternOrigin=pattern.copy()
            self.cellParameter=np.array([self.data['_cell_length_a'],self.data['_cell_length_b'],self.data['_cell_length_c'],self.data['_cell_angle_alpha'],self.data['_cell_angle_beta'],self.data['_cell_angle_gamma']],dtype=float)         
            self._out2cif_(position,label,'ligand')
            print('ligand structure generated')
            index_H_bonded=index_H
            for i in range(0,np.shape(sym)[0]):
                xT=eval(sym.iloc[i]['x'])
                yT=eval(sym.iloc[i]['y'])
                zT=eval(sym.iloc[i]['z'])    
                temp=np.stack((xT,yT,zT), axis=-1)
                position=np.append(position[:],temp,axis=0)
                label=np.append(label[:],labelOrigin.copy(),axis=0)
                atomType=np.append(atomType[:],atomTypeOrigin.copy(),axis=0)
                pattern=np.append(pattern[:],patternOrigin.copy(),axis=0)
                index_H_bonded=np.append(index_H_bonded,index_H_bonded,axis=0)
           # for i in range(0,len(index_H_bonded)):
           #     if index_H_bonded[i]:
           #         index_H_bonded[i]+=i
                    
            #self.position=np.expand_dims(position, axis=0)
            self._uniquePos_(position)
            positionFromCif=np.delete(position,self.dup_label,0)[:,:]
            
            labelFromCif=np.delete(label,self.dup_label,0)
            atomTypeFromCif=np.delete(atomType,self.dup_label,0)
            patternFromCif=np.delete(pattern,self.dup_label,0)
            self.index_H_bonded=index_H_bonded
            
            self._out2cif_(positionFromCif,atomTypeFromCif,'Rigid')
            print('rigid structure generated')
            
            return positionFromCif,labelFromCif,atomTypeFromCif,patternFromCif
        
        def _Occupancy_(self,inputProperty):
            deletList=np.where(self.occupancy<0.5)
            outputProperty=np.delete(inputProperty,[deletList],0)
            return outputProperty
        
        def _uniquePos_(self,position):
            self.dup_label=[]
            positionPBC=self._pbc_(position)
            for i in range(0,np.shape(position)[0]):
                for j in range (0,i):
                    if np.array_equal(positionPBC[i,:],positionPBC[j,:]):
                        self.dup_label.append(i)
            self.dup_label=np.unique(self.dup_label)
            
            return self.dup_label
        
        
        def _out2cif_(self,position=None,atomNames=None,outputFilename=None,cellParameter=None):
            if atomNames is None:
                atomNames=self.atomNames
            if cellParameter is None:
                cellParameter=self.cellParameter
            var=DataOut()
            var._out2cif_('cifBasedOut/'+self.cifName,atomNames, cellParameter, position, 'none',outputFilename)
            
            return
        
        def _cellGen_(self):
            
            #for exp_cif
            M=np.zeros([6])
            M[0:3]=self.cellParameter[0:3]
            M[3:6]=self.cellParameter[3:6]*np.pi/180

            omega=M[0]*M[1]*M[2]*np.sqrt(1-np.cos(M[3])**2-np.cos(M[4])**2-np.cos(M[5])**2+2*np.cos(M[3])*np.cos(M[4])*np.cos(M[5]))
            self.transM=[[M[0],M[1]*np.cos(M[5]),M[2]*np.cos(M[4])],
                     [0,M[1]*np.sin(M[5]),M[2]*(np.cos(M[3])-np.cos(M[4])*np.cos(M[5]))/np.sin(M[5])],
                     [0,0,omega/(M[0]*M[1]*np.sin(M[5]))]]

            # for clean
            M=np.zeros([6])
            M[0:3]=np.array(self.structure_clean.lattice.abc,dtype=float)
            M[3:6]=np.array(self.structure_clean.lattice.angles,dtype=float)*np.pi/180
            self.cellParameterCoRE=M.copy()
            self.cellParameterCoRE[3:6]=self.cellParameterCoRE[3:6]*180/np.pi
            omega=M[0]*M[1]*M[2]*np.sqrt(1-np.cos(M[3])**2-np.cos(M[4])**2-np.cos(M[5])**2+2*np.cos(M[3])*np.cos(M[4])*np.cos(M[5]))
            self.transMCoRE=[[M[0],M[1]*np.cos(M[5]),M[2]*np.cos(M[4])],
                     [0,M[1]*np.sin(M[5]),M[2]*(np.cos(M[3])-np.cos(M[4])*np.cos(M[5]))/np.sin(M[5])],
                     [0,0,omega/(M[0]*M[1]*np.sin(M[5]))]]
            
            return
        
        def _frac2Cartesian_(self,positionIn):
            self._cellGen_()
            positionOut=np.matmul(self.transM,(positionIn).T).T
            
            return positionOut
        
        def _frac2CartesianCoRE_(self,positionIn):
            self._cellGen_()
            positionOut=np.matmul(self.transMCoRE,(positionIn).T).T
            
            return positionOut
        
  
        def _fullRank_add_Gaussian_(self,positionR,pattern,index_associatedH,index_H,transM,imgOutNum):
            temp=np.arange(np.shape(positionR)[0])
            index_other=np.setdiff1d(temp,np.array(index_H))
            modifiedPositionFract=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            modifiedPositionCartn=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            for i in range(0,imgOutNum):
                for k in range(0,np.shape(index_other)[0]):
                    vib=np.random.normal(0,pattern[index_other[k]])*self.T_ratio
                    modifiedPositionCartn[i,index_other[k]]=positionR[index_other[k]]+vib
                    if index_associatedH[index_other[k]] < np.shape(positionR)[0]:
                           modifiedPositionCartn[i,index_associatedH[index_other[k]]]=positionR[index_associatedH[index_other[k]]]+vib        
                modifiedPositionFract[i]=np.matmul(np.linalg.inv(transM),modifiedPositionCartn[i].T).T
            return modifiedPositionFract,modifiedPositionCartn
        
        def _fullRank_add_Gaussian_aniso(self,positionR,index_associatedH,index_H,transM,imgOutNum):
            def ellipseGen(arrayIn,factor):
                arrayT=1.5874*np.array([[arrayIn[0],arrayIn[3],arrayIn[4]],[arrayIn[3],arrayIn[1],arrayIn[5]],[arrayIn[4],arrayIn[5],arrayIn[2]]])
                array=np.linalg.inv(arrayT)
                flag=False        
                while flag==False:
                    x=float(2*np.sqrt(arrayIn[0])*(np.random.rand(1,1)-0.5))
                    y=float(2*np.sqrt(arrayIn[1])*(np.random.rand(1,1)-0.5))
                    z=float(2*np.sqrt(arrayIn[2])*(np.random.rand(1,1)-0.5))               
                    v=np.array([x,y,z])
                    if np.matmul(np.matmul(v.T,array),v) < factor**2:
                        flag=True          
                return v
            def ellipseGen_multiGaussian(arrayIn,factor):
                cov=np.array([[arrayIn[0],arrayIn[3],arrayIn[4]],[arrayIn[3],arrayIn[1],arrayIn[5]],[arrayIn[4],arrayIn[5],arrayIn[2]]])
                v=np.random.multivariate_normal([0,0,0],cov)
                return v*factor
            
            temp=np.arange(np.shape(positionR)[0])
            index_other=np.setdiff1d(temp,np.array(index_H))
            modifiedPositionFract=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            modifiedPositionCartn=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            for i in range(0,imgOutNum):
                for k in range(0,np.shape(index_other)[0]):
                    patternArray=self.pattern_aniso[np.where(self.pattern_aniso_label==self.labelFromCif[index_other[k]])[0][0],:]
                    vib=ellipseGen_multiGaussian(patternArray,self.T_ratio)
                    modifiedPositionCartn[i,index_other[k]]=positionR[index_other[k]]+vib
                    if index_associatedH[index_other[k]] < np.shape(positionR)[0]:
                           modifiedPositionCartn[i,index_associatedH[index_other[k]]]=positionR[index_associatedH[index_other[k]]]+vib        
                modifiedPositionFract[i]=np.matmul(np.linalg.inv(transM),modifiedPositionCartn[i].T).T
            return modifiedPositionFract,modifiedPositionCartn
        
        def _findHList_(self,pos,atomType):
            index_H=np.where(atomType=='H')
            temp=np.arange(np.shape(pos)[0])
            index_other=np.setdiff1d(temp,np.array(index_H[0]))
            varH=scipy.spatial.KDTree(pos[index_H])
            varOther=scipy.spatial.KDTree(pos[index_other])
            result=varH.sparse_distance_matrix(varOther,2)
            index_nonzero=result.nonzero()
            
            #########based on H, but output in heavy atom#########
            index_associatedH=(np.shape(pos)[0])*np.ones(np.shape(pos)[0]+1,dtype=int)
            distanceArray=2*np.ones(len(index_H[0]),dtype=float)
            index_HAssociated2=(np.shape(pos)[0])*np.ones(len(index_H[0]),dtype=int)
            if len(index_nonzero[0]) < len(index_H[0]):
                print('may have serious H neignbour error!!')
            for i in range(0,len(index_nonzero[0])):
                if result[index_nonzero[0][i],index_nonzero[1][i]] < distanceArray[index_nonzero[0][i]]:
                    distanceArray[index_nonzero[0][i]]=result[index_nonzero[0][i],index_nonzero[1][i]]                 
                    index_HAssociated2[index_nonzero[0][i]]=index_other[index_nonzero[1][i]]
            for i in range(0,len(index_H[0])):
                index_associatedH[index_HAssociated2[i]]=index_H[0][i]
            
            return index_associatedH,index_H[0]

        def _pbc_(self,pos):
            posNew=pos-0.5
            posOut=posNew-np.round(posNew)+0.5
            return posOut
        
        def _quickSuperCell_(self,pos):
            for i in range(0,3):
                pos_1=pos.copy()
                pos_1[:,i]=pos[:,i]+1
                pos_2=pos.copy()
                pos_2[:,i]=pos[:,i]-1
                pos=np.append(pos,pos_1,axis=0)
                pos=np.append(pos,pos_2,axis=0)
            return pos
      
        def _matchStructure_(self,positionFromCif,labelFromCif,atomTypeFromCif,structure_clean):
             """
             match Metal
             """
             metalAtomIndex=np.where(np.array(structure_clean.species,dtype=str)=='Zn')
             position_Mclean=structure_clean.frac_coords[metalAtomIndex]
             position_Mclean=self._frac2CartesianCoRE_(self._pbc_(position_Mclean))
             
             metalAtomIndex2=np.where(atomTypeFromCif=='Zn')
             position_M=self._frac2Cartesian_(self._quickSuperCell_(self._pbc_(positionFromCif[metalAtomIndex2])))
             

             #d=distance.cdist(np.mean(position_M,axis=0).reshape(1,3),position_M,'euclidean')
             #centerMIndex=np.argsort(d)[0,0:np.shape(position_Mclean)[0]]
             #displacement=np.mean(position_Mclean,axis=0)-np.mean(position_M[centerMIndex],axis=0)
             #position_Mclean=position_Mclean-displacement
             #d=distance.cdist(np.mean(position_Mclean,axis=0).reshape(1,3),position_Mclean,'euclidean')
             #position_McleanSorted=np.array([ position_Mclean[x] for i,x in enumerate(np.argsort(d))][0])
             ### y: cif position ; x: CoRE position
             
             def max_cost1(w,x=position_Mclean,y=position_M):
                def nearest_neighbors_kd_tree(x, y, k) :
                    x, y = map(np.asarray, (x, y))
                    tree =scipy.spatial.cKDTree(y)
                    ordered_neighbors = tree.query(x, k)[1]
                    nearest_neighbor = np.empty((len(x),), dtype=np.intp)
                    nearest_neighbor.fill(-1)
                    used_y = set()
                    for j, neigh_j in enumerate(ordered_neighbors) :
                        for k in neigh_j :
                            if k not in used_y :
                                nearest_neighbor[j] = k
                                used_y.add(k)
                                break
                    return nearest_neighbor

                r=R.from_euler('xyz',w[0:3])      
                xNew = r.apply(x)+w[3:6]
                
                neignbour=nearest_neighbors_kd_tree(xNew, y, 2)
                distArray=np.array([xNew[i]-y[neignbour[i]] for i in range(0,np.shape(x)[0])],dtype=float)
                dist=np.sum(distArray**2)
                print(dist)
                return dist
             
             def recursion(w,i):
                 result = minimize(max_cost1, w).fun  
                 i=i+1
                 if result < 0.5 and i <2 :
                     r=R.from_euler('xyz',result.x[0:3])
                     position_Mclean_final=r.apply(position_Mclean)+result.x[3:6]
                 else:
                     wNew=[np.random.random_sample((3,))*360,np.random.random_sample((3,))*5]
                     recursion(wNew,i)
                     print('+++++++++++'+str(i)+'++++++++++++')
                     
                 return position_Mclean_final
             
             a=recursion(np.array([1,1,1,0,0,0]),0)

             
             #result = minimize(max_cost2, w)
             #r=R.from_euler('xyz',result.x)
             
            
             return a
           
             
#             position_clean=structure_clean.cart_coords[:]
#             matchedStructure=position_clean
#             
#             return matchedStructure
#       
#
#             def max_cost2(w,x=position_McleanSorted,y=np.array(position_M[centerMIndex])):
#                r=R.from_euler('xyz',w)      
#                xNew = r.apply(x)
#                dist=0
#                distNow=0
#                tuple1 = tuple(list(range(0,np.shape(y)[0]))) 
#                for i in range(0,np.shape(y)[0]):
#                    distNow=distance.cdist(xNew[i].reshape(1,3),[y[x] for i,x in enumerate(tuple1)],'euclidean')**2
#                    distMin=np.amin(distNow)
#                    tupleT=list(tuple1)
#                    tupleT.remove(tupleT[int(np.where(distNow==distMin)[1])])
#                    tuple1=tuple(tupleT)
#                    dist=dist+distMin
#                    
#                print(dist)
#                return dist
        
 