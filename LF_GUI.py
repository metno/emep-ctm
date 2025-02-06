#!/usr/bin/python3

# Written by Johan S. Wind April 2017
# Modified by Peter Wind June 2017

# For ubuntu install libraries with:
# sudo apt-get install python3-pyqt5 
# sudo apt install python3-netcdf4 

# requires "world_50m.txt" which can be googled

import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from math import *
import re
import netCDF4

class ContourPlotter:
    def __init__(self, par):
        self.par = par
        f = open("world_50m.txt")
        t = f.read()
        f.close()
        self.data = []
        proj = self.par.fh.projection
        print('projection: ',proj)
        for i in t.split("\n\n"):#split into blocks of "islands"
            self.data.append([])#list of lists
            for j in i.split("\n"):
                if not j: continue
                x, y = j.split()
                if(proj == 'Stereographic'):
                    xp = 8*50
                    yp = 110*50
                    fi = -32.0
                    an = 237.73*50
                    if float(y)<0: y='0'
                    x_ps = an*tan(pi*0.25-radians(float(y))*0.5)*sin(radians(float(x)-fi))
                    y_ps = -an*tan(pi*0.25-radians(float(y))*0.5)*cos(radians(float(x)-fi))
                    self.data[-1].append([x_ps, y_ps])
                else:
                    self.data[-1].append([float(x), float(y)])
        self.convert2pixels()

    def convert2pixels(self):
        latname = self.par.fh.variables[self.par.reg_poll].dimensions[2]
        lonname = self.par.fh.variables[self.par.reg_poll].dimensions[3]
        lat = self.par.fh.variables[latname][:]
        lon = self.par.fh.variables[lonname][:]

        lata = 1./(lat[1]-lat[0])
        latb = 0.5-lata*lat[0]

        lona = 1./(lon[1]-lon[0])
        lonb = 0.5-lona*lon[0]
        
        self.lines = []
        for i in self.data:
            self.lines.append([])
            for j in i:
                x, y = j
                lati = lata*y+latb
                loni = lona*x+lonb
                if loni >= 0 and lati >= 0 and loni < len(lon) and lati < len(lat):
                    self.lines[-1].append([loni, lati])
                else:
                    self.lines.append([])

class DimEdit(QWidget):
    def __init__(self, par, dimid, map_or_over, AreaActived):
        super(DimEdit, self).__init__(par)#init QWidget

        self.par = par #parent
        self.name = self.par.fh.variables[self.par.reg_poll].dimensions[dimid] # fh file handle
        vlayout = QVBoxLayout() #init boxes over each others
        vlayout.addWidget(QLabel(self.name)) #add a text to widget
        mbutton = QPushButton("-") #init a button
        mbutton.setMaximumSize(30, 30) 
        mbutton.clicked.connect(self.decrement) #self.decrement is called when mbutton clicked
        pbutton = QPushButton("+")
        pbutton.setMaximumSize(30, 30)
        pbutton.clicked.connect(self.increment)
        self.valueLabel = QLabel()
        adjust = QHBoxLayout()#init boxes on side of eachothers
        adjust.addWidget(mbutton)
        adjust.addWidget(self.valueLabel)
        adjust.addWidget(pbutton)
        if self.name == "klevel": 
            self.par.dim[map_or_over][dimid] = self.par.fh[self.par.reg_poll].shape[dimid]-1


        vlayout.addLayout(adjust)#put the QHBox in the QVBox
        self.setLayout(vlayout)#activate

        self.dimid = dimid
        self.setMinimumSize(180, 60)
        self.setMaximumSize(180, 60)
        self.map_or_over = map_or_over
        self.AreaActived = AreaActived
        self.updateText()

    def updateText(self):#write text in side bar
        self.valueLabel.setText("%d"%self.par.dim[self.map_or_over][self.dimid])
        try: #if string value is defined
            self.valueLabel.setText(str(int(self.par.fh[self.name][self.par.dim[self.map_or_over][self.dimid]])))
        except:
            pass

    def decrement(self):
        self.par.dim[self.map_or_over][self.dimid] = max(self.par.dim[self.map_or_over][self.dimid]-1, 0)#in side bar
        if self.map_or_over == 0:
            self.par.mapView.updateImage()#in map
        else:
            self.par.mapView.updateOverlay()#in map
        self.updateText()

    def increment(self):
        maxsize = self.par.fh[self.par.reg_poll].shape[self.dimid]-1
        self.par.dim[self.map_or_over][self.dimid] = min(self.par.dim[self.map_or_over][self.dimid]+1, maxsize)#in side bar
        if self.map_or_over == 0:
            self.par.mapView.updateImage()
        else:
            self.par.mapView.updateOverlay()
        self.updateText()

def func(v):
#    return v**.3 #NB: must change ifunc consisently
    return v**.2
def ifunc(v):
    return v**5
    
def toColor(v):
    x=1
    if(isnan(270-func(v)*270)):
        x=0
    else:
        x=int(270-func(v)*270)
    return QColor.fromHsv(x, 255, 255)#from [0 1] to rainbow
def toColor_nofunc(v): #don't scale v first
    return QColor.fromHsv(int(270-v*270), 255, 255)

def myfloat2str(v):#nicer format
    if v < 1e-2:
        return ("%.1e"%v).replace("e-0", "e-")
    else:
        return "%.3f"%v

class ColorScale(QWidget): #scale bar
    
    def __init__(self, par, scale=1.0):
        super(ColorScale, self).__init__(par)
        self.par = par
        self.setMinimumHeight(50)
        self.setMaximumHeight(50)
        self.repaint()
        self.scale=scale
        self.unit=""
        self.name=""

    def paintEvent(self, e):
        p = QPainter()
        p.begin(self)
        w, h = self.frameGeometry().width()-20, self.frameGeometry().height()
        p.setPen(QColor(0,0,0))#black
        for i in range(w):
            p.fillRect(i, 0, 1, 30, toColor_nofunc((i+.5)/w))#make rainbow slicewise
        if(self.scale < 1.02 and self.scale>0.98) : self.scale=1
        for i in range(10): #put numbers
            v = (i+.5)/10
            if i==9: v=1.0
            x = int(w*v)
            p.drawLine(x, 20, x, 30)#tick
            p.drawText(QRectF(x-50, 30, 100, 20), Qt.AlignCenter, myfloat2str(self.scale*ifunc(v)))
        p.drawText(QRectF(0.85*w, 5, 100, 20), Qt.AlignCenter, self.unit)
        p.drawText(QRectF(0.5*w, 0, 200, 20), Qt.AlignCenter, self.name)
        p.end()

class MapView(QWidget):
    def __init__(self, par):
        super(MapView, self).__init__(par)
        self.par = par
        self.data = []

        s_r = self.par.fh.variables[self.par.reg_poll].shape#reg all dimensions 
        s_l = self.par.fh.variables[self.par.loc_poll].shape#loc all dimensions 
        self.data_h, self.data_w = s_r[2], s_r[3]
        if len(s_l)>5:
            self.overh, self.overw = s_l[4], s_l[5]

        self.img = QImage(self.data_w, self.data_h, QImage.Format_ARGB32)#empty image container

        self.updateImage()

#        self.resize(self.data_w*4, self.data_h*4)

        self.overlay = 0
        self.ActivatedM = [[0 for col in range(self.overh)] for row in range(self.overw)]
        #self.overx = self.overy = 0

        self.contour = ContourPlotter(self.par)#coastal line

    def updateImage(self):
        if  len(self.par.fh.variables[self.par.reg_poll].dimensions)==6:
            self.data = self.par.fh.variables[self.par.reg_poll][self.par.dim[0][0],self.par.dim[0][1],:,:,self.overw//2,self.overh//2]
        elif  len(self.par.fh.variables[self.par.reg_poll].dimensions)==4:
            self.data = self.par.fh.variables[self.par.reg_poll][self.par.dim[0][0],self.par.dim[0][1],:,:]
        else:
            print("ERROR")
    
        maxval = max(self.data.flatten())

        #if re.search(r'fraction', self.par.reg_poll): 
        #    oscale = 1.0
        #else:
        #    oscale = 1.0/(1.E-6+ maxval)
        oscale = 1.0/(1.E-6+ maxval)
            
        for y in range(self.data_h):
            for x in range(self.data_w):
                v = self.data[y,x]
                v *= oscale

                self.img.setPixel(x, self.data_h-1-y, toColor(v).rgb())#put data in img as color
        self.par.colorScale_reg.scale = maxval
        self.unit_reg = self.par.fh.variables[self.par.reg_poll].units
        self.par.colorScale_reg.unit = self.par.fh.variables[self.par.reg_poll].units
        self.par.colorScale_reg.name = self.par.fh.variables[self.par.reg_poll].long_name
        self.par.colorScale_reg.update()
        self.repaint()

    def updateOverlay(self):
        if self.par.AreaActivated and self.overlay:
            adata = self.par.fh.variables[self.par.loc_poll][self.par.dim[1][0],self.par.dim[1][1],self.data_h-1-self.par.areay,self.par.areax,:,:]
            dx = self.par.areax - self.par.overx
            dy = self.par.areay - self.par.overy
            active_sign = 1
            if (abs(dx)<self.overw/2 and abs(dy)<self.overh/2):
                if self.ActivatedM[int(dy+self.overh/2)][int(dx+self.overw/2)] == 0:
                    self.ActivatedM[int(dy+self.overh/2)][int(dx+self.overw/2)] = 1
                    for y in range(self.overh):
                        if y+dy < self.overh:
                            for x in range(self.overw):
                                if x+dx < self.overw:
                                    self.par.odata[y,x] += adata[y+dy,x+dx]
                                
                else:
                    self.ActivatedM[int(dy+self.overh/2)][int(dx+self.overw/2)] = 0
                    for y in range(self.overh):
                        if y+dy < self.overh:
                            for x in range(self.overw):
                                if x+dx < self.overw:
                                    self.par.odata[y,x] -= adata[y+dy,x+dx]
                    
        else:
            self.par.odata = self.par.fh.variables[self.par.loc_poll][self.par.dim[1][0],self.par.dim[1][1],self.data_h-1-self.par.overy,self.par.overx,:,:]
            self.ActivatedM = [[0 for col in range(self.overh)] for row in range(self.overw)]           
            self.ActivatedM[int(self.overh/2)][int(self.overw/2)] = 1

        print(self.par.fh.variables[self.par.loc_poll].long_name)

        self.overh, self.overw = self.par.odata.shape
        self.overimg = QImage(self.overw, self.overh, QImage.Format_ARGB32)

        self.res = 1
        numbers = re.findall("[^0-9][0-9]+x([0-9]+)", "@"+self.par.fh.variables[self.par.loc_poll].long_name+"@")
        if len(numbers):
            self.res = int(numbers[-1])
        print("source size:",self.res)

#        self.res = 1
#        numbers = re.findall("[^0-9]([0-9]+)[^0-9]", "@"+self.par.reg_poll+"@")
#        if len(numbers):
#            self.res = int(numbers[-1])
#        print(numbers)

        sum = 0
        oscale = 1.0
        Nsum = 0
        for y in range(self.overh):
            for x in range(self.overw):
                Nsum += self.ActivatedM[y][x]

        maxval = max(self.par.odata.flatten())
        if re.search(r'transport', self.par.loc_poll): oscale = 1.0/(1.E-6+self.par.odata[int((self.overh-1)/2),int((self.overw-1)/2)])

        oscale = 1.0/(1.E-6+ maxval)

        for y in range(self.overh):
            for x in range(self.overw):
                v = self.par.odata[self.overh-1-y,x]
                sum += v
                v *= oscale
                self.overimg.setPixel(x, y, toColor(v).rgb())
        self.par.slabel.setText("Local sum: %.1f%%"%((sum*100)/Nsum))
        self.par.colorScale_loc.scale = 1.0/oscale
        self.par.colorScale_loc.unit = self.par.fh.variables[self.par.loc_poll].units
        self.par.colorScale_loc.name = self.par.fh.variables[self.par.loc_poll].long_name

        self.par.colorScale_loc.update()
        self.repaint()

    def mouseMoveEvent(self, event):
        x, y = event.pos().x(), event.pos().y()
        self.overx, self.overy = x*self.data_w//self.frameGeometry().width(), y*self.data_h//self.frameGeometry().height()
        proj = self.par.fh.projection
        if(proj == 'lon lat'):
            latname = self.par.fh.variables[self.par.reg_poll].dimensions[2]
            lonname = self.par.fh.variables[self.par.reg_poll].dimensions[3]
            lat = self.par.fh.variables[latname][int(self.data_h-1-self.overy)]
            lon = self.par.fh.variables[lonname][int(self.overx)]
        else:
            lat = self.par.fh.variables['lat'][int(self.data_h-1-self.overy),int(self.overx)]
            lon = self.par.fh.variables['lon'][int(self.data_h-1-self.overy),int(self.overx)]            
        txt = 'lat='+str(lat)+' lon='+str(lon)
        self.setToolTip(txt)
    def mousePressEvent(self, e):
        if self.overlay:
            if self.par.AreaActivated:
                x, y = e.pos().x(), e.pos().y()
                self.par.areax, self.par.areay = x*self.data_w//self.frameGeometry().width(), y*self.data_h//self.frameGeometry().height()
                self.setMouseTracking(1)
                self.updateOverlay()
            else:
                self.overlay = 0
                self.ActivatedM = [[0 for col in range(self.overh)] for row in range(self.overw)]
                self.setMouseTracking(1)
                self.repaint()
        else:
            x, y = e.pos().x(), e.pos().y()
            self.par.overx, self.par.overy = x*self.data_w//self.frameGeometry().width(), y*self.data_h//self.frameGeometry().height()
            self.par.slabel.setText("Local sum:")
            self.overlay = 1
            self.setToolTip('')
            self.setMouseTracking(0)
            self.ActivatedM[int(self.overh/2)][int(self.overw/2)] = 1
            self.updateOverlay()

    def paintEvent(self, e):
        p = QPainter()
        p.begin(self)
        p.drawImage(QRect(0, 0, self.frameGeometry().width(), self.frameGeometry().height()), self.img)

        w, h = self.frameGeometry().width(), self.frameGeometry().height()
        sw = w/self.data_w
        sh = h/self.data_h
        pen = QPen(QColor(0,0,0))#for coastal line
        pen.setWidth(2)
        p.setPen(pen)
        for i in self.contour.lines:
            last = 0
            for j in i:
                if last:
                    p.drawLine(int(last[0]*sw), int(h-last[1]*sh), int(j[0]*sw),int( h-j[1]*sh))
                last = j

        if self.overlay or self.par.AreaActivated:            
            rect = QRectF(self.par.overx*sw-(self.overw-1)/2*sw*self.res, self.par.overy*sh-(self.overh-1)/2*sh*self.res, self.overw*sw*self.res, self.overh*sh*self.res)
            rect2 = QRect(int(rect.left()-3), int(rect.top()-3), int(rect.width()+7), int(rect.height()+7))

            p.fillRect(rect2, QColor(0,0,0))
            p.drawImage(rect, self.overimg)
            for y in range(self.overh):
                for x in range(self.overw):
                    if(self.ActivatedM[y][x]):
                        p.drawText(int((self.par.overx+x-self.overw/2+1)*sw-w/2), int((self.par.overy+y-self.overh/2+1)*sh-h/2), int(w),int(h), Qt.AlignCenter, 'x')

        p.end()

class Main(QWidget):
    def __init__(self, fn):
        super(Main, self).__init__()
        vlayout = QVBoxLayout()
        hlayout = QHBoxLayout()
        sidebar = QVBoxLayout()

        self.dim = [[0 for i in range(6)] for j in range(2)] #define a 6x2 matrix

        self.colorScale_reg = ColorScale(self, 1.0)
        self.colorScale_loc = ColorScale(self, 1.0)

        self.filename = fn
        try:
            self.fh = netCDF4.Dataset(self.filename, mode="r")
        except:
            print ("Error opening file")
            self.openFile()

        for var in self.fh.variables:
            if re.search(r'fraction', var): 
                self.reg_poll = var #init
                self.loc_poll = var #init
                break

        abutton = QPushButton("AreaSum")
        self.AreaActivated = 0
        abutton.clicked.connect(self.AreaActivate)
        self.areax = 0
        self.areay = 0

        self.slabel = QLabel("Local sum: ")
        font = QFont()
        font.setPixelSize(20)
        font.setBold(True)
        self.slabel.setFont(font)

        obutton = QPushButton("Open file")
        obutton.clicked.connect(self.openFile)

        qbutton = QPushButton("Quit")
        qbutton.clicked.connect(self.close)

        self.slabel.setMaximumSize(180, 30)
        obutton.setMaximumSize(180, 30)
        qbutton.setMaximumSize(180, 30)
        abutton.setMaximumSize(180, 30)

        self.adjusts = [DimEdit(self, [0,1][i%2], i//2, 0) for i in range(4)]

        reg_spec = QComboBox(self)
        print(len(self.fh.variables))
        for var in self.fh.variables:
            if len(self.fh.variables[var].dimensions)==6: reg_spec.addItem(var)
            if len(self.fh.variables[var].dimensions)==4: reg_spec.addItem(var)

        reg_spec.setMaximumSize(180, 30)

        sidebar.addWidget(QLabel("Regional"))

        reg_spec.activated[str].connect(self.set_pollutant_reg)

        sidebar.addWidget(reg_spec)

        for i in range(2):
            sidebar.addWidget(self.adjusts[i])


        sidebar.addWidget(QLabel("Urban"))
        sidebar.addWidget(abutton)
 
        sidebar.addWidget(self.slabel)

        loc_spec = QComboBox(self)
        for var in self.fh.variables:
            if re.search(r'fraction', var): loc_spec.addItem(var)

        loc_spec.activated[str].connect(self.set_pollutant_loc)

        sidebar.addWidget(loc_spec)

        for i in range(2):
            sidebar.addWidget(self.adjusts[i+2])

        self.mapView = MapView(self)

        sidebar.addWidget(obutton)
        sidebar.addWidget(qbutton)

        vlayout.addWidget(self.colorScale_reg)
        vlayout.addWidget(self.colorScale_loc)
        vlayout.addWidget(self.mapView)

        hlayout.addLayout(vlayout)
        hlayout.addLayout(sidebar)
        self.setLayout(hlayout)
        
        self.resize(900, 900)
        self.move(1100,200)
        self.show()


    def set_pollutant_reg(self, text):
        self.reg_poll = text
        self.mapView.updateImage()
    def set_pollutant_loc(self, text):
        self.loc_poll = text
        self.mapView.updateOverlay()

    def AreaActivate(self):
        if(self.AreaActivated):
            self.AreaActivated = 0
            print("deactivated")
        else:
            self.AreaActivated = 1
            print("activated")
        

    def openFile(self):
        fname, idontknow = QFileDialog.getOpenFileName(self, "Open file", ".")
        if not fname: return

        self.filename = fname
        try:
            #We could save this if it is a performance issue
            self.fh = netCDF4.Dataset(self.filename, mode="r")
            self.close()
            Main(self.filename)
        except:
            print ("Error opening file")
            self.openFile()

def main():
    app = QApplication([])
    app.setStyle("cleanlooks")
    name = "uEMEP_full.nc"
    if len(sys.argv) <= 1:
        print ("No file given, defaulting to", name)
    else:
        name = sys.argv[1]
    win = Main(name)
    sys.exit(app.exec_())

if __name__ == "__main__": main()
