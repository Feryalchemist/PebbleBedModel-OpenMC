# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 07:49:30 2021

@author: ferya
"""

import math
import os
import openmc
import openmc.deplete
from random import seed, random, choice
from pyevtk.hl import pointsToVTK
from threading import Thread

def vtk_PolyDataFile(data,directory):
    import numpy as np
    for key in data:
        data[key]=np.array(data[key])
    pointsToVTK(directory,
                data['x'],data['y'],data['z'],
                data=data)

# You can change the pebble grouping method
def radial_axial(Pebbledf,count,**kwargs):
    rd_size=kwargs['radial_segment']
    ax_size=kwargs['axial_segment']
    an_size=kwargs['angle_segment']
    Pebbledf = Pebbledf.copy()

    z = [Pebbledf['z'][i] for i in range(len(Pebbledf['id'])) if Pebbledf['type'][i]==1]
    p = [min(z)+(i*(max(z)-min(z))/ax_size) for i in range(0,ax_size+1)]
    a = [(i*(max(Pebbledf['rho'])-min(Pebbledf['rho']))/rd_size) for i in range(0,rd_size+1)]
    t = [(i*(2*math.pi)/an_size) for i in range(0,an_size+1)]
    
    label = [0]*len(Pebbledf['id'])
    repeat=True
    count+=1
    
    bound=[]
    for i in range(0,rd_size*ax_size*an_size):
        bound.append([int(i/an_size/rd_size)%ax_size,
                      int(i/an_size)%rd_size,
                      i%an_size,
                      count+i])
    
    for k in range(0,len(Pebbledf['type'])):
        for q in range(0,len(bound)):
            if (Pebbledf['type'][k]!=2 and
                Pebbledf['z'][k] > p[bound[q][0]] and 
                Pebbledf['z'][k] <= p[bound[q][0]+1] and
                Pebbledf['rho'][k] > a[bound[q][1]] and
                Pebbledf['rho'][k] <= a[bound[q][1]+1] and
                Pebbledf['theta'][k] > t[bound[q][2]] and
                Pebbledf['theta'][k] <= t[bound[q][2]+1]):
                label[k]=bound[q][3]                
    return(label)

def individual_grouping(Pebbledf,count,**kwargs):
    Pebbledf=Pebbledf.copy()
    label = [0]*len(Pebbledf['id'])
    count+=1
    for i in range(0,len(Pebbledf['id'])):
        if Pebbledf['type'][i]!=2:
            label[i]=count
            count+=1
    return(label)

# You can use different variation of TRISO distribution 
def sparse_TRISO(n=8335,r_TRISO=0.04075,r_PEBBLE=2.5,r_BOUND=0.1):    
    seed(123457)
    cnt=0
    rd=[[],[],[]]
    r=r_PEBBLE-(r_TRISO)
    rm=r_BOUND
    center=[0,0,0]
    trial=0

    def sphere(array,center):
        rr=math.sqrt(((array[0]-center[0])**2)+((array[1]-center[1])**2)+((array[2]-center[2])**2))
        return rr

    while (cnt<8335):
        x=float((random()*r*2)-r)
        y=float((random()*r*2)-r)
        z=float((random()*r*2)-r)
        rad=sphere([x,y,z],center)**2 - r**2
        test=0
        if (rad < 0):
            if cnt>0:
                for i in range (0,len(rd[0])):
                    rds = (sphere([x,y,z],[rd[0][i],rd[1][i],rd[2][i]])**2) - rm**2
                    if rds > 0:
                        test+=1
                        if test==len(rd[0]):
                            rd[0].append(x); rd[1].append(y); rd[2].append(z)             
                            cnt+=1
                trial+=1
            else:
                cnt+=1
                rd[0].append(x); rd[1].append(y); rd[2].append(z)
                
    return {'x':rd[0],'y':rd[1],'z':rd[2]}    

def random_ordered_TRISO(r_TRISO=0.04075, center=[0]*3, pitch=0.15, size = 27, combination = 25, randomize = True):
    if size%2!=1:
        size+=1
    if isinstance(combination,str):
        combination = size**3
    c_vec=[[],[],[]]
    if randomize==True:
        for i in range(0,combination):
            c_vec[0]+=random()*(pitch-r_TRISO)*2-(pitch-r_TRISO)
            c_vec[1]+=random()*(pitch-r_TRISO)*2-(pitch-r_TRISO)
            c_vec[2]+=random()*(pitch-r_TRISO)*2-(pitch-r_TRISO)
    else:
        c_vec[0]=[0]*combination; c_vec[1]=[0]*combination; c_vec[2]=[0]*combination
    pos={'x':[],'y':[],'z':[]}
    for i in range(0,size):
        for j in range(0,size):
            for k in range(0,size):
                pos['x']+=[center[0]-int(size/2)*pitch+i*pitch]+c_vec[0][int(random()*size)]
                pos['y']+=[center[1]-int(size/2)*pitch+j*pitch]+c_vec[1][int(random()*size)]
                pos['z']+=[center[2]-int(size/2)*pitch+k*pitch]+c_vec[2][int(random()*size)]
    return pos

class DEM_data():
    def __init__(self,step,total_step):
        self.step = step
        self.__store = [0]*total_step
        
    def label(self,function,**kwargs):
        if self.step==0:
            self.Pebbledf['label'] = function(Pebbledf=self.Pebbledf,count=0,**kwargs)
        else:
            self.Pebbledf['label'] = self.__store.copy()[self.step-1]['label']
            unlabelled = {}
            for key in self.Pebbledf:
                if key!='label':
                    unlabelled[key] = self.Pebbledf.copy()[key][len(self.Pebbledf['label']):len(self.Pebbledf['id'])]
            self.Pebbledf['label'] += function(Pebbledf = unlabelled,
                                               count = max(self.Pebbledf['label']),**kwargs)
    
    def read(self,filepath,unit_factor=100):
        point_position = endpoint_position = id_position = endid_position=2e+20
        type_position = endtype_position = vel_position = endvel_position=1e+20
        omega_position = endomega_position = radii_position = endradii_position=1e+20
        
        radiiList=[]
        pointList=[]
        idList=[]
        typeList=[]
        velList=[]
        omegaList=[]
        
        def polar_angle(x,y):
            theta=math.atan(y/x)
            if x>=0 and y>=0:
                theta=theta
            if x<0 and y>=0:
                theta=theta+math.pi
            if x<0 and y<0:
                theta=theta+math.pi
            if x>=0 and y<0:
                theta=theta+math.pi*2
            return theta 
    
        with open(filepath) as fp:
            cnt=0
            for line in fp:
                for part in line.split():
                    if "POINTS" in part:
                        point_data=line.replace('POINTS ','').replace(' float\n','')
                        point_data=int(point_data)
                        point_position=cnt
                        endpoint_position=math.ceil(point_position+point_data/3)
                    elif "id" in part:
                        id_data=line.replace('id 1 ','').replace(' int\n','')
                        id_data=int(id_data)
                        id_position=cnt
                        endid_position=math.ceil(id_position+id_data/9)
                    elif "type" in part:
                        type_data=line.replace('type 1 ','').replace(' int\n','')
                        type_data=int(type_data)
                        type_position=cnt
                        endtype_position=math.ceil(type_position+type_data/9)
                    '''
                    elif "radius" in part:
                        radii_data=line.replace('radius 1 ','').replace(' double\n','')
                        radii_data=int(radii_data)
                        radii_position=cnt
                        endradii_position=math.ceil(radii_position+radii_data/9)
                    elif "omega" in part:
                        omega_data=line.replace('omega 3 ','').replace(' double\n','')
                        omega_data=int(omega_data)
                        omega_position=cnt
                        endomega_position=math.ceil(omega_position+omega_data/3)
                    elif "v 3" in line:
                        vel_data=line.replace('v 3 ','').replace(' double\n','')
                        vel_data=int(vel_data)
                        vel_position=cnt
                        endvel_position=math.ceil(vel_position+vel_data/3)
                    '''
                    
                
                if (cnt>point_position and cnt<endpoint_position+1):
                    pointList.extend(line.replace(' \n','').split(' '))
                elif (cnt>id_position and cnt<endid_position+1):
                    idList.extend(line.replace(' \n','').split(' '))
                elif (cnt>type_position and cnt<endtype_position+1):
                    typeList.extend(line.replace(' \n','').split(' '))
                '''
                elif (cnt>radii_position and cnt<endradii_position+1):
                    radiiList.extend(line.replace(' \n','').split(' '))
                elif (cnt>vel_position and cnt<endvel_position+1):
                    velList.extend(line.replace(' \n','').split(' '))
                elif (cnt>omega_position and cnt<endomega_position+1):
                    omegaList.extend(line.replace(' \n','').split(' '))
                '''
                cnt+=1
        
        ids=[int(i) for i in idList]
        self.Pebbledf = {'x':[],'y':[],'z':[],'vx':[],'vy':[],'vz':[],'wx':[],'wy':[],'wz':[],
                         'id':[],'rho':[],'theta':[],'r':[],'type':[]}
        for i in range(1,len(ids)+1):
            j = ids.index(i)
            self.Pebbledf['x'].append(float(pointList[j*3+0])*unit_factor)
            self.Pebbledf['y'].append(float(pointList[j*3+1])*unit_factor)               
            self.Pebbledf['z'].append(float(pointList[j*3+2])*unit_factor)
            '''
            self.Pebbledf['vx'].append(float(velList[j*3+0])*unit_factor)
            self.Pebbledf['vy'].append(float(velList[j*3+1])*unit_factor)
            self.Pebbledf['vz'].append(float(velList[j*3+2])*unit_factor)
            self.Pebbledf['wx'].append(float(omegaList[j*3+0]))
            self.Pebbledf['wy'].append(float(omegaList[j*3+1]))
            self.Pebbledf['wz'].append(float(omegaList[j*3+2]))
            self.Pebbledf['r'].append(float(radiiList[j])*unit_factor)
            '''
            self.Pebbledf['id'].append(int(idList[j]))
            self.Pebbledf['rho'].append(math.sqrt((float(pointList[j*3+0])**2)+
                                                  (float(pointList[j*3+1])**2)))
            self.Pebbledf['theta'].append(polar_angle(float(pointList[j*3+0]),
                                                      float(pointList[j*3+1])))
            self.Pebbledf['type'].append(int(typeList[j]))
        self.__store[self.step] = self.Pebbledf

class Recirculation():
    def __init__(self,step, timeline, burn_power, chain):
        self.step = step
        self.chain = openmc.deplete.Chain.from_xml(chain)
        self.timeline = timeline
        self.power = burn_power
        self.__storeBOC = [0]*(len(timeline)+1)
        self.__storeEOC = [0]*len(timeline)
        self.prefix_material_number = 300
        
    def read(self,filepath):
        results = openmc.deplete.ResultsList.from_hdf5(filepath)
        depleted_mat = results[-1].data.shape[1]
        self.matlist = [0]*depleted_mat

        summ=[]; suma=[]
        for i in range(0,depleted_mat):
            self.matlist[i] = {'atom den':[], 'atom fr':[],'density':0,'volume':0,
                          'mass':[],'mass fr':[], 'material_group':0, 'isotope':[]}
            summ +=[sum([results[-1].data[0,i,results[-1].nuc_to_ind[k]]*
                         int(''.join(filter(str.isdigit,str(k)))[:3])/
                         (6.023*10e23) for k in results[-1].nuc_to_ind])]
            suma +=[sum(results[-1].data[0,i,:])]

            self.matlist[i]['volume'] = results[-1].volume[str(i+1)]
            self.matlist[i]['material_group'] = i+1
            self.matlist[i]['density'] = summ[i]/results[-1].volume[str(i+1)]

            for n in results[-1].nuc_to_ind:
                resa = results[-1].data[0,i,results[-1].nuc_to_ind[n]]
                mrel = ''.join(filter(str.isdigit,str(n)))[:3]
        
                self.matlist[i]['isotope'].append(str(n))
                self.matlist[i]['atom den'].append(resa*1e-24/self.matlist[i]['volume'])
                self.matlist[i]['atom fr'].append(resa/suma[i])
                self.matlist[i]['mass'].append(resa*int(mrel)/(6.023*10e23))
                self.matlist[i]['mass fr'].append(resa*int(mrel)/(6.023*10e23)/summ[i])

        self.__storeBOC[self.step+1] = self.matlist
        self.__storeEOC[self.step]   = self.matlist 
    
    def add(self,isotope,fraction,density,volume,material_group,Pebbledf,fraction_type='atom fr'):
        material_framelist = ({'isotope':isotope,'mass':[0]*len(isotope),'activity':[0]*len(isotope),
                               'spec act':[0]*len(isotope),'atom den':[0]*len(isotope),'atom fr':[],'mass fr':[],'material_group':material_group,'density':density,'volume':Pebbledf['label'].count(material_group)*volume})
        #Please review this part
        if fraction_type == 'mass fr':
            material_framelist['mass fr']=fraction
            sum_atom_fraction = sum([fraction[i]/(int(''.join(filter(str.isdigit,str(isotope[i])))[:3])) for i in range(0,len(isotope))])
            material_framelist['atom fr']=[fraction[i]/(int(''.join(filter(str.isdigit,str(isotope[i])))[:3]))/sum_atom_fraction for i in range(0,len(isotope))]
            
        if fraction_type == 'atom fr':
            material_framelist['atom fr']=fraction
            sum_mass_fraction = sum([fraction[i]*(int(''.join(filter(str.isdigit,str(isotope[i])))[:3])) for i in range(0,len(isotope))])
            material_framelist['mass fr']=[fraction[i]*(int(''.join(filter(str.isdigit,str(isotope[i])))[:3]))/sum_mass_fraction for i in range(0,len(isotope))]
            
        if density>=0:
            material_framelist['atom den']=[density*fraction[i] for i in range(0,len(isotope))]

        else:
            material_framelist['atom den']=[density*fraction[i]*0.6023/(int(''.join(filter(str.isdigit,str(isotope[i])))[:3])) for i in range(0,len(isotope))]

        if self.__storeBOC[self.step]==0:
            self.__storeBOC[self.step] = [material_framelist]
            
        elif self.__storeBOC[self.step]!=0:
            locator = [x['material_group'] for x in self.__storeBOC[self.step]]
            if material_framelist['material_group'] in locator:
                a = self.__storeBOC[self.step][locator.index(material_framelist['material_group'])]
                b = material_framelist
                
                similar_za = list(set(b['isotope']).intersection(set(a['isotope'])))
                differ_za  = list(set(b['isotope'])-set(a['isotope']))
                for K in similar_za:
                    for key in a.keys():
                        if key != 'material_group' and key != '' and key != 'density':
                            a[key][a['isotope'].index(K)] += b[key][b['isotope'].index(K)]
                for K in differ_za:
                    for key in a.keys():
                        if key != 'material_group' and key != 'density':
                            a[key].append(b[key][b['isotope'].index(K)])
                a['mass fr'] = [a['mass fr'][x]/sum(a['mass fr']) for x in a['mass fr']]
                a['atom fr'] = [a['atom fr'][x]/sum(a['atom fr']) for x in a['atom fr']]
                a['density'] = sum(a['atom den'])
                self.__storeBOC[self.step][locator.index(material_framelist['material_group'])] = a
            else:
                self.__storeBOC[self.step].append(material_framelist)
        
        self.matlist = self.__storeBOC[self.step]
            
    def write_material(self,fraction_type='atom fr'):
        self.kernel_material = []
        if fraction_type == 'atom fr':
            tog = ['ao','atom/b-cm']
        else:
            tog = ['wo','g/cm3']
        for num,mat in enumerate(self.matlist,1):
            self.kernel_material += [openmc.Material(name='Kernel_'+str(int(num)))]
            self.kernel_material[num-1].set_density(tog[1],mat['density'])
            self.kernel_material[num-1].volume = mat['volume']
            for mum,iso in enumerate(mat['isotope'],0):
                self.kernel_material[num-1].add_nuclide(iso,mat[fraction_type][mum],tog[0])       
        
        
    
