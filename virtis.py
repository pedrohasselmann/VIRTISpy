#!/usr/bin/python
# -*- coding: utf-8 -*-

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

__version__="1.0.0"

__doc__='''
  ############################ VIRTIS READER ##############################
  ### Filename: virtis.py
  ### Author: VIRTIS-VEX / Romolo Politi (INAF-Italy)
  ### Adapted and debugged by Pedro H. Hasselmann for Rosetta/VIRTIS-M-VIS-IR/H
  ###
  ###
  ### Load VIRTIS spectro-image cubes
  ###
  #########################################################################
''' 
#import numba as nb
'''spec = [
('header',dict),
('dataType',str),
('cube',nb.float32[:]),
('table',nb.float32[:])
]

@nb.jitclass(spec)'''
class read_VIRTIS(object):
    """
        PDS data reader for VIRTIS (Visible and Infrared Thermal Imaging Spectrometer) on board of the mission Venus Express
    """
    from numpy import float32, transpose, array
    import struct, gc
    global gc, array, struct, transpose, float32

    def __init__(self,fn='',mission='ROSETTA', getHeader=True, getdata=True): # VEX or ROSETTA
        self.mission = mission
        self.fn = fn
        self.readHeader()
        if getdata is True: self.readPDSFile()
        if getHeader is False: self.header = None

    def readHeader(self):
       """
         Internal function. Read the input file and extract the PDS header
       """
       self.header={}
       in_file = open(self.fn,"rb")
       in_file.seek(0,2)
       
       if in_file.tell() == 0: # Empty file
         return
       else:
         in_file.seek(0)

       go_on = 1
       while 1:
         line = in_file.readline().split()
         #print(line)
         if len(line) == 1 and line[0] == 'END': # End of header
              break
         if not line: # Empty line
              go_on = 1
              continue
         if line[0] == '/*': # Comment
              continue
         if line[0] == 'NOTE':
              go_on = 0
              continue
         if line[0][0] == '"' or line[0][0].isdigit():
                for s in line: self.header[key] += s
                continue                      
         if go_on == 1:
                self.header[line[0]] = ' '.join(line[2:])
                key = line[0]

       in_file.close()
       #raw_input(self.header.keys())
       line, key = None, None

    def getBand(self,band):
       if self.dataType == 'VIRTIS DATA':
              ret=self.cube[band,:,:]
       elif self.dataType== 'VIRTIS GEOMETRY':
              ret=self.cube[band,:,:]*self.geoCoeff[band]
       return ret.astype(float32)

    def getValue(self,band,sample,line):
       if self.dataType == 'VIRTIS DATA':
              ret=self.cube[band,sample,line]
       elif self.dataType== 'VIRTIS GEOMETRY':
              ret=self.cube[band,sample,line]*self.geoCoeff[band]
       return ret

    def getSpectrum(self,sample,line):
      if self.dataName=='VIRTIS_H':
         return self.cube[:,sample]
      else:
         return self.cube[:,sample,line]

    def getWave(self):
       if self.dataName=='VIRTIS_H':
         return self.table[0,:]
       else:
         return self.table[:,0,0]

    def getSCET(self,line):
       return array(self.HK[0:3,0,line])

    def getHK(self,val,line):
       if type(val).__name__=='str':
         if val in self._HKName:
              idx=self.HKName.index(val)
         else:
              raise Exception(val+' is not a valid Housekeeping record')

       else:
         if val < len(self.HKName):
              idx=val
         else:
              raise Exception('Not valid Housekeeping record >0 <%d'%(len(self._HKName)))

       if type(line).__name__== 'list':
         ret=[]
         for i in line:
              ret.append(self.HK[idx,0,i])
       else:
         ret=self.HK[idx,0, line]
       return ret

    def readPDSFile(self):
       """
          Data reader
          Supported Formats:
         -- VIRTIS M VIS
         -- VIRTIS M IR Raw
         -- VIRTIS M IR Calibrated
         -- VIRTIS M IR Geometry
         -- VIRTIS H
       """
       header = self.header
       recByte=int(header['RECORD_BYTES'])
       dataPos=int(header['^QUBE'])-1  # Data Pointer
       dim=map(int,header['CORE_ITEMS'][1:-1].split(','))
       sfxItem=map(int,header['SUFFIX_ITEMS'][1:-1].split(','))
       crType=header['CORE_ITEM_TYPE']

       if crType == 'REAL' or crType == '"REAL"':
         tp='f'
         crByte=4
       elif crType == 'MSB_INTEGER':
         crByte=int(header['CORE_ITEM_BYTES'])
         if crByte ==4:
              tp='i'
         elif crByte == 2:
              tp='h'
       
       
       self.dataName = header[self.mission+':CHANNEL_ID'][1:-1]
       self.dataType=header['STANDARD_DATA_PRODUCT_ID'][1:-1]
       procLev=int(header['PROCESSING_LEVEL_ID'])
       #print(self.dataName, self.dataType)
       
       if self.dataType=='VIRTIS DATA' and self.dataName!='VIRTIS_H':

         coreBlock=(dim[0]+sfxItem[0])*(dim[1]+sfxItem[1])*dim[2]*crByte
         shape=(dim[0]+sfxItem[0],(dim[1]+sfxItem[1]),dim[2])
         lineDim=dim[2]*crByte

         dataStruc=str(dim[0])+tp+sfxItem[0]*'L'
         d=(dim[1]+sfxItem[1])*dim[2]*dataStruc

         fl=open(self.fn,'rb')
         fl.seek(dataPos*recByte)
         buff=fl.read()
         out=struct.unpack_from('>'+d,buff)
         
         #print(dataStruc)
         #print('dimension',dim,crType,shape[0]*shape[2])
         
         if procLev==3:
              #scetDim=2
              temp=[out[i] for i  in range(dim[0], len(out), dim[0]+1)]
              temp=array(filter(lambda a: a !=0,temp), dtype=float32)

              fl.seek((dataPos*recByte)+coreBlock)
              buff=fl.read()
              bott=struct.unpack_from(('>'+(3*dim[0]*dim[1])*'f'),buff)
              
              self.cube=array(out, dtype=float32).reshape(shape,order='F')[0:dim[0],:,:]
              self.scet=temp[1::3]
              self.table=array(bott, dtype=float32).reshape((dim[0],dim[1],3),order='F')
              bott, buff, temp = None, None, None

         elif procLev==2:
              data=array(out, dtype=float32).reshape(shape,order='F')
              self.cube=data[:,0:-1,:]
              tmp=list(data[:,-1,:].reshape((dim[0]*dim[2]),order='F'))
              tmp2=struct.pack(len(tmp)*'h',*tmp)
              hk=struct.unpack_from(len(tmp)*'H',tmp2)
              tmp=array(hk).reshape((dim[0],dim[2]),order='F')
              tmp2=tmp[0:(82*5),:]
              self.HK=tmp2.reshape((82,5,dim[2]),order='F')
              self.HKName=['Data SCET-1','Data SCET-2','Data SCET-3','Acquisition ID','# of subslices + 1st serial #','Data Type','SPARE','ME_default HK SCET-1','ME_default HK SCET-2','ME_default HK SCET-3','V_MODE','ME_PWR_STAT','ME_PS_TEMP','ME_DPU_TEMP','ME_DHSU_VOLT','ME_DHSU_CURR', 'EEPROM_VOLT','IF_ELECTR_VOLT','SPARE','M_ME_general HK SCET-1','M_ME_general HK SCET-2','M_ME_general HK SCET-3','M_ECA_STAT','M_COOL_STAT', 'M_COOL_TIP_TEMP','M_COOL_MOT_VOLT','M_COOL_MOT_CURR','M_CCE_SEC_VOLT','SPARE','MVIS_HK_report SCET-1','MVIS_HK_report SCET-2','MVIS_HK_report SCET-3','M_CCD_VDR_HK','M_CCD_VDD_HK','M_+5_VOLT','M_+12_VOLT','M_-12_VOLT','M_+20_VOLT','M_+21_VOLT','M_CCD_LAMP_VOLT','M_CCD_TEMP_OFFSET','M_CCD_TEMP','M_CCD_TEMP_RES','M_RADIATOR_TEMP','M_LEDGE_TEMP','OM_BASE_TEMP','H_COOLER_TEMP','M_COOLER_TEMP','M_CCD_WIN_X1','M_CCD_WIN_Y1', 'M_CCD_WIN_X2','M_CCD_WIN_Y2','M_CCD_DELAY','M_CCD_EXPO','M_MIRROR_SIN_HK','M_MIRROR_COS_HK','M_VIS_FLAG_ST','SPARE','MIR_HK_report SCET-1', 'MIR_HK_report SCET-2','MIR_HK_report SCET-3','M_IR_VDETCOM_HK','M_IR_VDETADJ_HK','M_IR_VPOS','M_IR_VDP','M_IR_TEMP_OFFSET','M_IR_TEMP', 'M_IR_TEMP_RES','M_SHUTTER_TEMP','M_GRATING_TEMP','M_SPECT_TEMP','M_TELE_TEMP','M_SU_MOTOR_TEMP','M_IR_LAMP_VOLT','M_SU_MOTOR_CURR', 'M_IR_WIN_Y1','M_IR_WIN_Y2','M_IR_DELAY','M_IR_EXPO','M_IR_LAMP_SHUTTER','M_IR_FLAG_ST','SPARE']
         else:
              print(procLev)
              raise Exception('No case')

       elif self.dataType=='VIRTIS GEOMETRY':

         fl=open(self.fn,'rb')
         fl.seek(dataPos*recByte)
         coreBlock=(dim[0]+sfxItem[0])*(dim[1]+sfxItem[1])*dim[2]*crByte
         buff=fl.read()

         shape=(dim[0]+sfxItem[0],(dim[1]+sfxItem[1]),dim[2])
         lineDim=dim[2]*crByte

         dataStruc=str(dim[0])+tp+sfxItem[0]*'L'
         d=(dim[1]+sfxItem[1])*dim[2]*dataStruc
         
         out=struct.unpack_from('>'+d,buff)
         
         #print('dimension',dim,sfxItem,crType,shape[0]*shape[2])

         #print('CORE_ITEMS: ',self._dim)
         geoCoeff=[0.0001]*100
         geoCoeff=array(geoCoeff)
         geoCoeff[13:14]=0.001
         geoCoeff[15]=0.00001
         geoCoeff[29]=0.001
         geoCoeff[32]=1
         #self._geoName=['Surf longit, corner1','Surf longit, corner2','Surf longit, corner3','Surf longit, corner4','Surf latit, corner1','Surf latit, corner2','Surf latit, corner3','Surf latit, corner4','Surf longit, center','Surf latit, center','Incidence at surf','Emergence at surf','Phase at surf','Elevation on surf layer','Slant distance','Local time','Cloud longit, corner1','Cloud longit, corner2','Cloud longit, corner3','Cloud longit, corner4','Cloud latit, corner1','Cloud latit, corner2','Cloud latit, corner3','Cloud latit, corner4','Cloud longit, center','Cloud latit, center','Incidence on clouds','Emergence on clouds','Phase on clouds','Elevation below clouds','Right ascension','Declination','M-common frame']
         self.geoCoeff=geoCoeff
         self.cube=array(out, dtype=float32).reshape(shape,order='F')


       elif self.dataName=='VIRTIS_H' and procLev ==3:

         tablePos=int(header['^TABLE'])-1 # Table Pointer

         coreBlock=(dim[0]+sfxItem[0])*(dim[1]+sfxItem[1])*dim[2]*crByte
         shape=(dim[0]+1,dim[2])#(dim[0]+sfxItem[0],dim[2])
         lineDim=dim[2]*crByte

         dataStruc=str(dim[0])+tp+'L'#sfxItem[0]*'L'
         d=(dim[1]+sfxItem[1])*(dim[2])*dataStruc
         
         # QUB Data
         fl=open(self.fn,'rb')
         fl.seek(dataPos*recByte)
         buff=fl.read()
         out=struct.unpack_from('>'+d,buff)
         self.cube=array(out, dtype=float32).reshape(shape,order='F')[0:dim[0],:]
         
         # Table - WAVE/FWHM/UNCERTAIN
         fl.seek(tablePos*recByte)
         buff=fl.read()
         bott=struct.unpack_from('>'+(3*dim[0])*'f',buff)
         self.table=array(bott, dtype=float32).reshape((3,dim[0]),order='F')

         buff, out, bott = None, None, None     

       fl.close()
       out = None
       gc.collect()

# END
