import pyglet as pg
from pyglet import shapes
from pyglet.gl import glClearColor
import math as mt
import numpy as np

def eng_string( x, format='%s', si=False):
    sign = ''
    if x < 0:
        x = -x
        sign = '-'
    exp = int( mt.floor( mt.log10( x)))
    exp3 = exp - ( exp % 3)
    x3 = x / ( 10 ** exp3)

    if si and exp3 >= -24 and exp3 <= 24 and exp3 != 0:
        exp3_text = 'yzafpnum kMGTPEZY'[ ( exp3 - (-24)) / 3]
    elif exp3 == 0:
        exp3_text = ''
    else:
        exp3_text = 'e%s' % exp3

    return ( '%s'+format+'%s') % ( sign, round(x3,1), exp3_text)

def getColor(colorVal, minColorVal, maxColorVal, colorFunc):
    try: x_ = float(colorVal - minColorVal)/(maxColorVal - minColorVal)
    except ZeroDivisionError: x_ = 0.5 # cmax == cmin

    x = colorFunc(x_)

    blue  = int(255* min((max((4*(0.75-x), 0.)), 1.)))
    red   = int(255* min((max((4*(x-0.25), 0.)), 1.)))
    green = int(255* min((max((4*mt.fabs(x-0.5)-1., 0.)), 1.)))
    return (red, green, blue)

class MeshRender:
    
    def __init__(self, Width = 1200, Height = 700, scale = 1):
        self.Width = Width
        self.Height = Height
        self.scale = scale
        self.center = True
        self.batch = pg.graphics.Batch()
        self.bkgColor = (255, 255, 255)
        self.barWidth = 1
        self.defaultElementColor = (255, 166, 234)
        self.barColor = (23, 23, 23)
        self.offsetx = 0
        self.offsety = 0
        self.legend = True
        self.legendDiscretize = 10
        self.nodesLabels = False
        self.labelAllNodes = False
        self.nodesSize = 2
        self.autoScale = False
        self.legendTitle = ''
        self.colorElements = True
        self.deform_scale = 1.0
    
    def drawElements(self, elements):
        
        self.window = pg.window.Window(self.Width,self.Height)
        glClearColor(self.bkgColor[0], self.bkgColor[1], self.bkgColor[2], 1.0)
        
        nodesGroup = pg.graphics.OrderedGroup(3)
        legendGroup = pg.graphics.OrderedGroup(2)
        barsGroup = pg.graphics.OrderedGroup(1)
        trianglesGroup = pg.graphics.OrderedGroup(0)
        
        bars = []
        triangles = []
        nodesPos = []
        nodes, nodesLabels = [],[]
        
        discret = self.legendDiscretize
        legenddx, legenddy = 20, 30
        legendx, legendy = self.Width - 5*legenddx, (self.Height - (legenddy*discret))/2
        legend = []
        
        if self.legend:
            ymin, ymax, xmin, xmax = elements[0].getminColorVal, elements[0].getmaxColorVal, 0, discret-1
            colorFunc = elements[0].getcolorFunc
            
            self.offsetx -= self.Width - 5*legenddx
            
            for i in range(discret):
                
                colorVal = (ymax- ymin)/(xmax-xmin)*float(i) + ymin
                
                y = legendy + i*legenddy
                legend.append(
                    shapes.Rectangle(x=legendx, y=y,width=legenddx, height=legenddy,
                                    color=getColor(colorVal, ymin, ymax, colorFunc),
                                    batch=self.batch,
                                    group=legendGroup
                                    )
                )
                
                legend.append(
                    pg.text.Label( '  ' + str(round(colorVal,0)),
                            font_name='Times New Roman',
                            font_size=10,
                            x=legendx + legenddx, y=y + legenddy/2,
                            color=(0, 0, 0, 255),
                            batch=self.batch,
                            group=legendGroup
                            )
                )
            
            legend.append(
                    pg.text.Label( self.legendTitle ,
                            font_name='Times New Roman',
                            font_size=10,
                            x=legendx, y= legendy + legenddy*(discret+0.5),
                            color=(0, 0, 0, 255),
                            batch=self.batch,
                            group=legendGroup
                            )
                )
            
        xmax, ymax, xmin, ymin = -9e9, -9e9, 9e9, 9e9
        for element in elements:
            for node in element.nodes:
                if node.x > xmax:
                    xmax = node.x
                if node.y > ymax:
                    ymax = node.y
                    
                if node.x < xmin:
                    xmin = node.x
                if node.y < ymin:
                    ymin = node.y
        
        if self.autoScale:
            
            minRatio = min([
                self.Width/abs(xmax - xmin),
                self.Height/abs(ymax - ymin)
            ])
            
            self.scale = minRatio*0.75
                        
        if self.center:
        
            self.offsetx = (self.Width - (xmax-xmin)*self.scale)/2 - xmin*self.scale
            self.offsety = (self.Height - (ymax-ymin)*self.scale)/2 - ymin*self.scale
        
        
        for element in elements:
            nodesPos.clear()
            for i in range(len(element.nodes)):
                
                j = (i+1)%3
                
                if self.deform_scale != 1.0:
                
                    x0 = self.offsetx + (element.nodes[i].x + element.nodes[i].dx*self.deform_scale)*self.scale 
                    y0 = self.offsety + (element.nodes[i].y + element.nodes[i].dy*self.deform_scale)*self.scale
                    x1 = self.offsetx + (element.nodes[j].x + element.nodes[j].dx*self.deform_scale)*self.scale
                    y1 = self.offsety + (element.nodes[j].y + element.nodes[j].dy*self.deform_scale)*self.scale
                
                else:
                    x0 = self.offsetx + (element.nodes[i].x)*self.scale 
                    y0 = self.offsety + (element.nodes[i].y)*self.scale
                    x1 = self.offsetx + (element.nodes[j].x)*self.scale
                    y1 = self.offsety + (element.nodes[j].y)*self.scale
                
                    
                l = ((x1 - x0)**2 + (y1 - y0)**2)**0.5
                angle = (180/mt.pi)*mt.asin((x1-x0)/l)

                bars.append(shapes.Line(x=x0, y=y0, x2=x1, y2=y1, width= self.barWidth, color=self.barColor,
                                    batch=self.batch, group=barsGroup))
                bars[-1].rotation = angle
                
                nodesPos.append([x0, y0])
                
                normalNode = False
                if self.nodesSize > 0:
                    if element.nodes[i].dfix:
                        nodeColor = (255, 48, 48)
                    elif element.nodes[i].externalForce:
                        nodeColor = (120, 255, 115)
                    else:
                        nodeColor = (150, 177, 255)
                        normalNode = True
                        
                    if not normalNode:
                        nodes.append(
                            shapes.Circle(x=x0, y=y0, radius=self.nodesSize, color=nodeColor,
                                        batch=self.batch, group=nodesGroup)
                        )
                        
                        nodes.append(
                            shapes.Arc(x=x0, y=y0, radius=self.nodesSize, color=(0, 0, 0), batch =self.batch, group=nodesGroup)
                        )
                        
                    if self.nodesLabels and (self.labelAllNodes or not normalNode):
                        nodesLabels.append(
                            pg.text.Label( str(element.nodes[i].id),
                                    font_name='Times New Roman',
                                    font_size=8,
                                    x=x0 + 5, y=y0 +5 ,
                                    color=(0, 0, 0, 255),
                                    batch=self.batch,
                                    group=nodesGroup
                                    )
                        )
                        
            if self.colorElements:
                color = element.getColor()
            else:
                color = self.defaultElementColor
                
            triangles.append(shapes.Triangle(
                x=nodesPos[0][0],
                y=nodesPos[0][1],
                x2=nodesPos[1][0],
                y2=nodesPos[1][1],
                x3=nodesPos[2][0],
                y3=nodesPos[2][1],
                color=color, batch=self.batch, group=trianglesGroup
                ))
            
                
        @self.window.event
        def on_draw():
            self.window.clear()
            self.batch.draw()
        
        pg.app.run()