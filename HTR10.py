#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 21:13:40 2021

@author: feryantama
"""
import math
import os
import openmc
import openmc.deplete
from random import seed, random, choice
from pyevtk.hl import pointsToVTK
from threading import Thread

class HTR10_Core_Universe():
    def __init__(self,Pebbledf,fuel_mat,**kwargs):
        #======================================================================
        self.Pebbledf=Pebbledf
        self.material=fuel_mat
        TRISOmatlist=[]
        mat =  {'nom':['buffer'   ,'IPYC'     ,'SiC'      ,'OPYC'     ,'Matrix'   ,'Dummy'   ],
                'den':[1.1        ,1.9        ,3.18       ,1.9        ,1.73       ,1.84      ],
                'C'  :[5.51511e-2 ,9.52610e-2 ,4.77597e-2 ,9.52610e-2 ,8.67377e-2 ,9.22528e-2],
                'B10':[1.58513e-8 ,2.73795e-8 ,0          ,2.73795e-8 ,2.49298e-8 ,2.54951e-9],
                'B11':[6.38035e-8 ,1.10206e-7 ,0          ,1.10206e-7 ,1.00345e-7 ,1.02621e-8],
                'Si' :[0          ,0          ,4.77597e-2 ,0          ,0          ,0         ]}
        for n in range(0,len(mat['nom'])):
            sums = sum([mat[key][n] for key in list(mat.keys())[2:]])
            mtr  = openmc.Material(name=mat['nom'][n])
            mtr.set_density('atom/b-cm',sums)
            mtr.add_element('C',  mat['C'][n]/sums)
            mtr.add_nuclide('B10',mat['B10'][n]/sums)
            mtr.add_nuclide('B11',mat['B11'][n]/sums)
            mtr.add_element('Si', mat['Si'][n]/sums)
            TRISOmatlist+=[mtr]
        #buffer_mat, IPYC_mat,SiC_mat,OPYC_mat,Matrix_mat,Dummy_mat = TRISOmatlist
        if 'coolant_mat' in kwargs:
            self.coolant_mat = kwargs['coolant_mat']
        else:
            self.coolant_mat = openmc.Material(name='Coolant')
            mats = [3.90968e-05,1.04885e-05+4.29014e-07,2.33838e-07,8.58028e-07]
            self.coolant_mat.set_density('atom/b-cm',sum(mats))
            self.coolant_mat.add_element('N' ,mats[0]/sum(mats),'ao')
            self.coolant_mat.add_element('O' ,mats[1]/sum(mats),'ao')
            self.coolant_mat.add_element('Ar',mats[2]/sum(mats),'ao')
            self.coolant_mat.add_element('H' ,mats[3]/sum(mats),'ao')
        #======================================================================
        if 'buffer_mat' in kwargs:
            self.buffer_mat = kwargs['buffer_mat']
        else:
            self.buffer_mat =TRISOmatlist[0]    
        if 'IPYC_mat' in kwargs:
            self.IPYC_mat = kwargs['IPYC_mat']
        else:
            self.IPYC_mat =TRISOmatlist[1]    
        if 'SiC_mat' in kwargs:
            self.SiC_mat = kwargs['SiC_mat']
        else:
            self.SiC_mat =TRISOmatlist[2]    
        if 'OPYC_mat' in kwargs:
            self.OPYC_mat = kwargs['OPYC_mat']
        else:
            self.OPYC_mat =TRISOmatlist[3]    
        if 'Matrix_mat' in kwargs:
            self.Matrix_mat = kwargs['Matrix_mat']
        else:
            self.Matrix_mat =TRISOmatlist[4]    
        if 'Dummy_mat' in kwargs:
            self.Dummy_mat = kwargs['Dummy_mat']
        else:
            self.Dummy_mat =TRISOmatlist[5]
        self.matr_collection=TRISOmatlist + fuel_mat +[self.coolant_mat]
     
    def create_TRISO(self):
        TRISOsurf=[]; TRISOcell=[]
        spherenum=[0.025,0.034,0.038,0.0415]
        TRISOsurf+=[openmc.Sphere(name='Kernel',r=spherenum[0])]
        self.TRISOuniverse=[]
        for n,mat in enumerate(self.material,1):
            TRISOcell+=[openmc.Cell(name='Kernel_'+str(n),fill=mat,region= -TRISOsurf[0])]
            for nom,m,mat in zip(['buffer','IPYC','SiC'],list(range(1,4)),
                                 [self.buffer_mat,self.IPYC_mat,self.SiC_mat]):
                TRISOsurf+=[openmc.Sphere(name=nom+'_'+str(n),r=spherenum[m])]
                TRISOcell+=[openmc.Cell(  name=nom+'_'+str(n),fill=mat, region= +TRISOsurf[m-1] & -TRISOsurf[m])]
            TRISOsurf+=[openmc.Sphere(name='OPYC_'+str(n),r=spherenum[-1])]
            TRISOcell+=[openmc.Cell(  name='OPYC_'+str(n),fill=self.OPYC_mat,region=+TRISOsurf[-1])]
            self.TRISOuniverse += [openmc.Universe(name='TRISO_'+str(n+1),
                                                  cells=TRISOcell[(n-1)*5:n*5])]
        self.surf_collection = TRISOsurf
        self.cell_collection = TRISOcell
        self.univ_collection = self.TRISOuniverse
        
    def create_Pebble(self,TRISO_distribution,rFZ=2.5,rTRISO=0.0455):
        self.rTRISO   = rTRISO #This default for HTR10
        self.rFZ      = rFZ    #This default for HTR10
        self.TRISOpos = TRISO_distribution
        self.Pebbleuniverse = []
        self.TRISOdf=TRISO_distribution
        
        for n,km in enumerate(self.material):
            iTRISOcell    = []
            iTRISOsurf    = []
            for m,i,j,k in zip(list(range(0,len(self.TRISOpos['x']))),
                               self.TRISOpos['x'],self.TRISOpos['y'],self.TRISOpos['z']):
                iTRISOsurf+=[openmc.Sphere(name='f'+str(n)+'T'+str(m),
                                           x0=i,y0=j,z0=k,r=rTRISO)]
                iTRISOcell+=[openmc.Cell(name='f'+str(n)+'T'+str(m),  
                                         fill=self.TRISOuniverse[n], region=-iTRISOsurf[m])]
                iTRISOcell[m].translation=(i,j,k)
                if m==0:
                    Matrix_region=+iTRISOsurf[m]
                else:
                    Matrix_region = Matrix_region & +iTRISOsurf[m]
            iTRISOsurf+=[openmc.Sphere(name='f'+str(n)+'Tmatrix_fill',r=rFZ)]
            iTRISOcell+=[openmc.Cell(name='f'+str(n)+'Tmatrix_fill', 
                                     fill=self.Matrix_mat, region= Matrix_region & -iTRISOsurf[-1])]
            iTRISOcell+=[openmc.Cell(name='f'+str(n)+'Tmatrix_shell', 
                                     fill=self.Matrix_mat, region= +iTRISOsurf[-1])]
            self.surf_collection+=iTRISOsurf
            self.cell_collection+=iTRISOcell
            self.Pebbleuniverse +=[openmc.Universe(name='Pebble_'+str(n),cells=iTRISOcell)]
            self.univ_collection+=[self.Pebbleuniverse[n]]
    
    def create_Core(self,rPebble=3):
        self.rPebble=rPebble
        cell = []; surf = [];
        for i,ids,typ,lbl,x,y,z in zip(list(range(0,len(self.Pebbledf['id']))),
                                       self.Pebbledf['id'],self.Pebbledf['type'],
                                       self.Pebbledf['label'],self.Pebbledf['x'],
                                       self.Pebbledf['y'],self.Pebbledf['z']):
            surf += [openmc.Sphere(name='P_'+str(ids),x0=x,y0=y,z0=z,r=rPebble)]
            if typ==1:
                cell += [openmc.Cell(name='P_'+str(ids)+'f',fill=self.Pebbleuniverse[lbl-1],region=-surf[i])]
                cell[i].translation=(x,y,z)
            else:
                cell += [openmc.Cell(name='P_'+str(ids)+'d',fill=self.Dummy_mat,region=-surf[i])]
            
            if i==0:
                Pebble_reg = +surf[i]
            else:
                Pebble_reg = Pebble_reg & +surf[i]
        cell += [openmc.Cell(name='coolant',fill =self.coolant_mat,region=Pebble_reg)]
        self.Coreuniverse = openmc.Universe(name='Core',cells=cell)
        self.cell_collection += cell
        self.surf_collection += surf
        self.univ_collection += [self.Coreuniverse]
    
    def visualize(self,r=90,h=[351.818-610,351.818-130],c=[0, 0, 46.818],
                  dirs=os.getcwd()+'/Template',res=1000):
        source = os.getcwd()
        if os.path.isdir(dirs)==False:
            os.makedirs(dirs)
        os.chdir(dirs)
        
        top = openmc.ZPlane(name='top',z0=h[1],boundary_type='vacuum')
        bot = openmc.ZPlane(name='bottom',z0=h[0],boundary_type='vacuum')
        cyl = openmc.ZCylinder(name='cylinder',x0=0,y0=0,r=r,boundary_type='vacuum')
        root_cell = openmc.Cell(name='root_cell',fill=self.Coreuniverse,region=-top & +bot & -cyl)
        root_univ = openmc.Universe(cells=[root_cell])
        
        geometry = openmc.Geometry()
        geometry.root_universe=root_univ
        geometry.export_to_xml()
        
        material = openmc.Materials(self.matr_collection)
        material.export_to_xml()
        
        settings = openmc.Settings()
        settings.batches = 550
        settings.inactive = 10
        settings.particles = int(100000)
        settings.material_cell_offsets=False
        settings.export_to_xml()

        plot1 = openmc.Plot(plot_id=1)
        plot1.filename = 'HTR10-core-xz'
        plot1.basis='xz'
        plot1.origin = c
        plot1.width = [r*2,h[1]-h[0]]
        plot1.pixels = [res,int(res*(h[1]-h[0])/(r*2))]
        plot1.color_by = 'material'

        plot2 = openmc.Plot(plot_id=2)
        plot2.filename = 'HTR10-core-xy'
        plot2.basis='xy'
        plot2.origin = c
        plot2.width = [r*2,r*2]
        plot2.pixels = [res,res]
        plot2.color_by = 'material'

        ty = choice(list(range(1,max(self.Pebbledf['label']))))
        fp = choice([i for i,t in zip(list(range(0,len(self.Pebbledf['id']))),self.Pebbledf['label']) if t == ty])
        x, y, z = [self.Pebbledf[l][fp] for l in ['x','y','z']]
        ft = choice(list(range(0,len(self.TRISOdf))))
        
        plot3 = openmc.Plot(plot_id=3)
        plot3.filename = 'Pebble'
        plot3.basis='xy'
        plot3.origin = [x,y,z]
        plot3.width = [self.rPebble*2,self.rPebble*2]
        plot3.pixels = [res,res]
        plot3.color_by = 'material'
        
        plot4 = openmc.Plot(plot_id=4)
        plot4.filename = 'TRISO'
        plot4.basis='xy'
        plot4.origin = [x+self.TRISOdf['x'][ft],y+self.TRISOdf['y'][ft],z+self.TRISOdf['z'][ft]]
        plot4.width = [self.rTRISO*2,self.rTRISO*2]
        plot4.pixels = [res,res]
        plot4.color_by = 'material'
        
        plots = openmc.Plots([plot1,plot2,plot3,plot4])
        plots.export_to_xml()
        
        openmc.plot_inline(plot1)
        openmc.plot_inline(plot2)
        openmc.plot_inline(plot3)
        openmc.plot_inline(plot4)
        
        #os.system('openmc --plot')
        os.chdir(source)
        
class HTR10_Control_Universe():
    def __init__(self):
        return
    def create_rod(self,step,**kwargs):
        if 'coolant_mat' in kwargs:
            coolant = kwargs['coolant_mat']
        else:
            coolant = openmc.Material(name='Coolant')
            mats = [3.90968e-05,1.04885e-05+4.29014e-07,2.33838e-07,8.58028e-07]
            coolant.set_density('atom/b-cm',sum(mats))
            coolant.add_element('N' ,mats[0]/sum(mats),'ao')
            coolant.add_element('O' ,mats[1]/sum(mats),'ao')
            coolant.add_element('Ar',mats[2]/sum(mats),'ao')
            coolant.add_element('H' ,mats[3]/sum(mats),'ao')

        Control_rod=openmc.Material(name='Control_rod')
        Control_rod.add_element('B',0.8,'ao')
        Control_rod.add_element('C',0.2,'ao')
        Control_rod.set_density('g/cm3', 1.7)
    
        Spacing=openmc.Material(name='Spacing')
        Spacing.add_element('Cr',0.18,'wo')
        Spacing.add_element('Fe',0.681,'wo')
        Spacing.add_element('Ni',0.1,'wo')
        Spacing.add_element('Si',0.01,'wo')
        Spacing.add_element('Mn',0.02,'wo')
        Spacing.add_element('C',0.001,'wo')
        Spacing.add_element('Ti',0.008,'wo')
        Spacing.set_density('g/cm3', 7.9)
    
        Caps=openmc.Material(name='Caps_n_Joints')
        Caps.add_element('Fe',1,'ao')
        Caps.set_density('atom/b-cm',0.04)        
        Z_control=[0,4.5,53.2,56.8,105.5,109.1,157.8,161.4,210.1,213.7,262.4,264.7]
        C_control=[2.75,2.95,3,5.25,5.3,5.5,6.499]
        M_control=[Spacing,coolant,
                   Control_rod,coolant,
                   Spacing,coolant]
        delt = (394.2-171.2)
        Z_base = 171.2+delt*step

        control_surf=list()
        contc=0
        for i in range(0,len(Z_control)):
            control_surf.append(openmc.ZPlane(z0=round(351.818+Z_control[-1]-(Z_control[i]+Z_base),3)))
            contc+=1

        for i in range(0,len(C_control)-1):
            control_surf.append(openmc.ZCylinder(x0=0.0,y0=0.0,r=C_control[i]))
            contc+=1

        control_cell=list()
        for i in range(0,5):
            for j in range(0,len(C_control)-1):
                if j!=len(C_control)-2:
                    reg_contr=+control_surf[len(Z_control)+j] & -control_surf[1+len(Z_control)+j] & -control_surf[1+i*2] & +control_surf[2+i*2]
                if j==len(C_control)-2:
                    reg_contr=+control_surf[len(Z_control)+j] & -control_surf[1+i*2] & +control_surf[2+i*2]
                control_cell.append(openmc.Cell(fill=M_control[j],region=reg_contr))
        for i in range(0,6):
            for j in range(0,2):
               if j==0:
                    control_cell.append(openmc.Cell(fill=Caps,region=(+control_surf[len(Z_control)] & 
                                                                      -control_surf[len(Z_control)+5] &
                                                                      -control_surf[i*2] & 
                                                                      +control_surf[1+i*2])))
               if j!=0:
                    control_cell.append(openmc.Cell(fill=coolant,region=(+control_surf[len(Z_control)+5] & 
                                                                                  -control_surf[i*2] & 
                                                                                  +control_surf[1+i*2])))
                  
        control_cell.append(openmc.Cell(fill=coolant,region=(-control_surf[len(Z_control)] &
                                                                      -control_surf[0] &
                                                                      +control_surf[len(Z_control)-1])))
        control_cell.append(openmc.Cell(fill=coolant,region=(+control_surf[0] | 
                                                                      -control_surf[len(Z_control)-1])))
        Control_univ=openmc.Universe()
        for cell in control_cell:
            Control_univ.add_cell(cell)
        self.control_matr = [Control_rod,Spacing,Caps,coolant]
        self.Z_base = Z_base
        self.control_cell = control_cell
        self.control_surf = control_surf
        self.control_univ = Control_univ
        
class HTR10_Reactor_Universe():
    def __init__(self,filepath,Pebbledf,CR_Universe,RC_Universe):
        self.filepath = filepath
        self.Pebbledf = Pebbledf
        self.CR_Universe = CR_Universe
        self.RC_Universe = RC_Universe
    def Read_Table_34(self):
        table_34 = {'Block':[],'atom den':[],'N':[],'O-part':[],'Ar':[],'H':[],'O-moist':[],'C':[],'B':[],'left':[],'right':[],'top':[],'bottom':[]}
        with open(self.filepath,'r') as fp:
            readline = fp.readlines()
            for num,line in enumerate(readline,0):
                if num!=0:
                    for m,key in enumerate(table_34.keys(),0):
                        if key != 'Block':
                            if line.split('\t')[m]!='-':
                                table_34[key] += [float(line.split('\t')[m])]
                            else:
                                table_34[key] += [line.split('\t')[m]]
                        else:
                            table_34[key] += [int(line.split('\t')[m])]
        table_34['O'] = [table_34['O-part'][i] + table_34['O-moist'][i] for i in range(len(table_34['Block']))]
        del table_34['O-part']
        del table_34['O-moist']
        self.table_34 = table_34
    def Create_Material(self):
        self.Reflector_mat = []
        for i in range (0,len(self.table_34['Block'])):
            self.Reflector_mat += [openmc.Material(name='Reflector-Block_'+str(int(i)))]
            self.Reflector_mat[i].set_density('atom/b-cm',self.table_34['atom den'][i])
            for atom in ['N','O','Ar','B','H','C']:
                if atom != 'Block' and atom != 'atom den':
                    self.Reflector_mat[i].add_element(atom,self.table_34[atom][i],'ao')
    def Create_Geometry(self):
        levellist=[]
        levelnum=[610.0,540.0,510.0,495.0,465.0,450.0,430.0,402.0,388.764,351.818,130.0,114.7,105.0,95.0,40.0,0]
        for i in range (0,len(levelnum)):
            level=openmc.ZPlane(name='level_'+str(i),z0=351.818-levelnum[i])
            levellist.append(level)
        radiallist=[]
        radialnum=[25.0,41.75,70.75,90.0,95.6,108.6,140.6,148.6,167.793,190.0]
        for i in range (0,len(radialnum)):
            radial=openmc.ZCylinder(name='cyl_'+str(i),x0=0.0,y0=0.0,r=radialnum[i])
            radiallist.append(radial)
        cone_center=-(90.0/((90.0-25.0)/(388.764-351.818)))
        cone=openmc.ZCone(x0=0.0,y0=0.0,z0=cone_center,r2=(90.0/-(cone_center))**2)
        if len(self.Pebbledf['r'])==0:
            rrr = 2.5
        else:
            rrr = sum(self.Pebbledf['r'])/len(self.Pebbledf['r'])
        LS=min(self.Pebbledf['z'])-rrr
        if LS<=-258.182:
            LS=-258.182
            
        levellist[0].boundary_type='vacuum'
        levellist[levelnum.index(0)].boundary_type='vacuum'
        radiallist[radialnum.index(190)].boundary_type='vacuum'

        low_surf=openmc.ZPlane(name='low_surf',z0=LS)
        low_surf7=openmc.ZPlane(name='cell7low',z0=LS-(540-495))
        low_surf81=openmc.ZPlane(name='cell81low',z0=LS-(610-495))
        self.Reflector_sur = [] + [levellist] + [radiallist] + [cone] + [low_surf,low_surf7,low_surf81]
        #============================== Create Cell ===============================
        Boring_reg = -radiallist[radialnum.index(190.0)]
        #///////////////////////////////////CORE////////////////////////////////////
        coreboundreg=((-radiallist[radialnum.index(90.0)] & -levellist[levelnum.index(130.0)] & +levellist[levelnum.index(351.818)]) | (-cone & -levellist[levelnum.index(351.818)] & +levellist[levelnum.index(388.764)]) |(-radiallist[radialnum.index(25.0)] & -levellist[levelnum.index(388.764)] & +low_surf ))
        coreboundcell=openmc.Cell(fill=self.RC_Universe,region=coreboundreg)
        #/////////////////////////////COOLANT CHANNEL//////////////////////////////
        Coolant_surf = []; Coolant_cell = []
        for i in range(0,20):
            deg = math.radians(9+i*18)
            Coolant_surf += [openmc.ZCylinder(name='Coolant_'+str(i),
                                              x0=144.6 * math.cos(deg),
                                              y0=144.6 * math.sin(deg),r=4)]
            Coolant_cell += [openmc.Cell(fill = self.Reflector_mat[5],region = -Coolant_surf[i] &
                                         -levellist[levelnum.index(105.0)] & 
                                         +levellist[levelnum.index(610.0)])]
            Boring_reg = Boring_reg & (+Coolant_surf[i] | 
                                       +levellist[levelnum.index(105.0)]|
                                       -levellist[levelnum.index(610.0)]) 
        #///////////////////////////////CONTROL ROD////////////////////////////////
        Control_surf = []; Control_cell = []
        for n,i in enumerate([k for k in range(0,20) if k not in [3,5,7,9,13,19,17]],0):
            deg = math.radians(9+i*18)
            Control_surf += [openmc.ZCylinder(name='Control_'+str(i),
                                              x0=102.1 * math.cos(deg),
                                              y0=102.1 * math.sin(deg),r=6.5)]
            C_Cell = openmc.Cell(fill = self.CR_Universe,region = -Control_surf[n] &
                                 -levellist[levelnum.index(0)] & 
                                 +levellist[levelnum.index(450.0)])
            C_Cell.translation=(102.1 * math.cos(deg),102.1 * math.sin(deg),0)
            Control_cell += [C_Cell]
            Boring_reg = Boring_reg & (+Control_surf[n] | 
                                       +levellist[levelnum.index(0)]|
                                       -levellist[levelnum.index(450.0)])
        #///////////////////////////ROUND KLAK CHANNEL/////////////////////////////
        RKLAK_surf = []; RKLAK_cell = []
        for n,i in enumerate([3,5,7,9,13,19,17],0):
            deg = math.radians(9+i*18)
            RKLAK_surf += [openmc.ZCylinder(name='RKLAK_'+str(i),
                                            x0=102.1 * math.cos(deg),
                                            y0=102.1 * math.sin(deg),r=3)]
            RKLAK_cell += [openmc.Cell(fill = self.Reflector_mat[5],region = -RKLAK_surf[n] &
                                       -levellist[levelnum.index(0)] & 
                                       +levellist[levelnum.index(130.0)])]
            RKLAK_cell += [openmc.Cell(fill = self.Reflector_mat[5],region = -RKLAK_surf[n] &
                                       -levellist[levelnum.index(388.764)] & 
                                       +levellist[levelnum.index(610.0)])]
            Boring_reg = Boring_reg & (+RKLAK_surf[n] | 
                                       +levellist[levelnum.index(0)]|
                                       -levellist[levelnum.index(130.0)])
            Boring_reg = Boring_reg & (+RKLAK_surf[n] | 
                                       +levellist[levelnum.index(388.764)]|
                                       -levellist[levelnum.index(610.0)])
        #////////////////////////////////OVAL KLAK/////////////////////////////////
        OKLAK_surf = []; OKLAK_cell = []
        klakcirc1=list()
        klakcirc2=list()
        klakplane1=list()
        klakplane2=list()
        klakwall1=list()
        klakwall2=list()
        
        for i,n in enumerate([3,5,7,9,13,17,19],0):
            dej = math.radians(9+n*18-90)
            deg = math.radians(9+n*18)
            dek = math.radians(9+n*18+90)

            a1,b1,c1,d1,n1 = [math.cos(deg),math.sin(deg),0,101.6,-1]
            a2,b2,c2,d2,n2 = [math.cos(deg),math.sin(deg),0,95.6,1]
            xr1,yr1        = [a1*98.6+math.cos(dek)*5,b1*98.6+math.sin(dek)*5]
            xr2,yr2        = [a1*98.6+math.cos(dej)*5,b1*98.6+math.sin(dej)*5]
            a3,b3,c3,d3,n3 = [math.cos(dek),math.sin(dek),0,5,0]
            a4,b4,c4,d4,n4 = [math.cos(dek),math.sin(dek),0,-5,0]
            
            KLAKdf = [deg,a1,b1,c1,d1,a2,b2,c2,d2,xr1,yr1,3,xr2,yr2,3,a3,b3,c3,d3,a4,b4,c4,d4]

            klakcirc1 += [openmc.ZCylinder(x0=KLAKdf[9],y0=KLAKdf[10],r=KLAKdf[11])]
            klakcirc2 += [openmc.ZCylinder(x0=KLAKdf[12],y0=KLAKdf[13],r=KLAKdf[14])]
            klakplane1+= [openmc.Plane(a=KLAKdf[1],b=KLAKdf[2],c=KLAKdf[3],d=KLAKdf[4])]
            klakplane2+= [openmc.Plane(a=KLAKdf[5],b=KLAKdf[6],c=KLAKdf[7],d=KLAKdf[8])] 
            klakwall1 += [openmc.Plane(a=KLAKdf[15],b=KLAKdf[16],c=KLAKdf[17],d=KLAKdf[18])]
            klakwall2 += [openmc.Plane(a=KLAKdf[19],b=KLAKdf[20],c=KLAKdf[21],d=KLAKdf[22])]

            OKLAK_cell += [openmc.Cell(fill=self.Reflector_mat[5],
                                      region=((-klakcirc1[i]| -klakcirc2[i]|
                                               (-klakplane1[i]& +klakplane2[i]& -klakwall1[i]& +klakwall2[i]))&
                                              (-levellist[levelnum.index(130.0)])&
                                              (+levellist[levelnum.index(388.764)])))]
            Boring_reg=Boring_reg & (((+klakcirc1[i]& +klakwall1[i])|
                                      (+klakcirc2[i]& -klakwall2[i])| +klakplane1[i]| -klakplane2[i])|
                                     +levellist[levelnum.index(130.0)]|
                                     -levellist[levelnum.index(388.764)])           
        OKLAK_surf += klakcirc1+klakcirc2+klakplane1+klakplane2+klakwall1+klakwall2
        #/////////////////////////////////HOT DUCT/////////////////////////////////
        Hduct_surf = [openmc.XCylinder(y0=0,z0=351.818-480,r=14.999),
                      openmc.XPlane(x0=90), openmc.XPlane(x0=190)]
        Hduct_cell = openmc.Cell(fill=self.Reflector_mat[5],name='Hotduct',region=(-Hduct_surf[0] & -Hduct_surf[2] & +Hduct_surf[1]))
        Boring_reg = Boring_reg & (+Hduct_surf[0] | +Hduct_surf[2] | -Hduct_surf[1])
        #//////////////////////////////REACTOR BLOCK///////////////////////////////
        CBlock_cell = []
        for i in range(0,len(self.table_34['Block'])):
            if self.table_34['Block'][i] == 83:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i]),fill=self.Reflector_mat[i],
                                            region=Boring_reg & -levellist[levelnum.index(351.818)] & -radiallist[radialnum.index(90.0)] & +cone & +levellist[levelnum.index(388.764)] & +radiallist[radialnum.index(25.0)])]
            if self.table_34['Block'][i] == 48:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i])+'a',fill=self.Reflector_mat[i],
                                            region=Boring_reg & -levellist[levelnum.index(40)] & -radiallist[radialnum.index(167.793)] & +levellist[levelnum.index(388.764)] & +radiallist[radialnum.index(148.6)])]
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i])+'b',fill=self.Reflector_mat[i],
                                            region=Boring_reg & -levellist[levelnum.index(40)] & -radiallist[radialnum.index(148.6)] & +levellist[levelnum.index(95)] & +radiallist[radialnum.index(108.6)])]
            if self.table_34['left'][i] == 0 and self.table_34['Block'][i]!=81 and self.table_34['Block'][i]!=7 and self.table_34['Block'][i]!=83 and self.table_34['Block'][i]!=48 and self.table_34['Block'][i]!=5:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i]),fill=self.Reflector_mat[i],
                                            region=Boring_reg & -levellist[levelnum.index(self.table_34['top'][i])] & -radiallist[radialnum.index(self.table_34['right'][i])] & +levellist[levelnum.index(self.table_34['bottom'][i])])]
            if self.table_34['left'][i] != 0 and self.table_34['Block'][i]!=81 and self.table_34['Block'][i]!=7 and self.table_34['Block'][i]!=83 and self.table_34['Block'][i]!=48 and self.table_34['Block'][i]!=5:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i]),fill=self.Reflector_mat[i],
                                            region=Boring_reg & -levellist[levelnum.index(self.table_34['top'][i])] & -radiallist[radialnum.index(self.table_34['right'][i])] & +levellist[levelnum.index(self.table_34['bottom'][i])] & +radiallist[radialnum.index(self.table_34['left'][i])])]
            if self.table_34['Block'][i]==81:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i]),fill=self.Reflector_mat[i],
                                            region=(-low_surf7 & -radiallist[radialnum.index(self.table_34['right'][i])] & +low_surf81))]
            if self.table_34['Block'][i]==7:
                CBlock_cell += [openmc.Cell(name='Block_'+str(self.table_34['Block'][i]),fill=self.Reflector_mat[i],
                                            region=(-low_surf  & -radiallist[radialnum.index(self.table_34['right'][i])] & +low_surf7))]
        self.Reflector_sur += Coolant_surf+Control_surf+RKLAK_surf+OKLAK_surf+Hduct_surf
        self.Reflector_cel =  CBlock_cell+Coolant_cell+Control_cell+RKLAK_cell+OKLAK_cell+[Hduct_cell]+[coreboundcell]
        
        Reactor_univ=openmc.Universe(cells=self.Reflector_cel)
        root_cell=openmc.Cell(name='root_cell',fill=Reactor_univ,region=-radiallist[radialnum.index(190)] & -levellist[levelnum.index(0.0)] & +levellist[levelnum.index(610.0)])
        self.root_universe=openmc.Universe(name='root_universe')
        self.root_universe.add_cell(root_cell)

