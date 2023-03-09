import numpy as np
import sys
from vispy import app,scene
from vispy.color import Color
from vispy.visuals.transforms import MatrixTransform
import time
from vispy.scene.visuals import Text
from vispy.visuals.filters import ShadingFilter
import math_dependencies as tr
phi=(1+np.sqrt(5))/2

canvas = scene.SceneCanvas(keys='interactive', size=(600, 600), show=True)

view = canvas.central_widget.add_view()
view.bgcolor = '#9f9f9f'
view.camera = 'turntable'
view.padding = 100

quad =np.array([
    [phi,phi-1,1],
    [1.5*phi-1,1.5*phi-1.5,phi-.5],
    [1,phi,phi-1],
    [phi-.5,1.5*phi-1,1.5*phi-1.5],
    ])/2

quads=[quad,quad@tr.get_transform([1,2,3,4,0]).T,quad@tr.get_transform([0,3,1,2,4]).T]

color=Color("#3f51b5")
shading_filter = ShadingFilter(shininess=100)

    
class TriangleBlock(scene.visuals.Mesh):
    def __init__(self,idcols,bkcol=(0,0,0),transform=np.eye(4)):
        V=np.array([
            [phi,1,phi-1],
            [1,phi-1,phi],
            [phi-1,phi,1],
            2*np.array([phi,1,phi-1])/3,
            2*np.array([1,phi-1,phi])/3,
            2*np.array([phi-1,phi,1])/3,
            [phi-.5,phi-.5,phi-.5]
        ])/2
        I=np.array([
            [0,1,6],
            [1,2,6],
            [2,0,6],
            [1,0,3],
            [1,3,4],
            [2,1,4],
            [2,4,5],
            [0,2,5],
            [0,5,4],
            [6,5,4],
    ])
        colors=list(idcols)+[(0,0,0)]*7
        scene.visuals.Mesh.__init__(self,vertices=V,faces=I,face_colors=colors,parent=view.scene,shading='flat')
        self.transform=MatrixTransform(transform)
        #self.direction=transform@np.array([1,1,1,1])


class PentagonBlock(scene.visuals.Mesh):
    def __init__(self,idcols,bkcol=(0,0,0),transform=np.eye(4)):
        v0=np.array([           
            [phi,1,phi-1],
            [1,phi-1,phi],
            [1,1-phi,phi],
            [phi,-1,phi-1],
            [2,0,0]])/2
        c=np.array([1,0,phi-1])
        V=np.vstack((v0,.5*v0+.5*c,2*v0/4))
        I=np.array([
            [0,1,5],
            [1,2,6],
            [2,3,7],
            [3,4,8],
            [4,0,9],
            [1,6,5],
            [2,7,6],
            [3,8,7],
            [4,9,8],
            [0,5,9],
            [1,0,10],
            [2,1,11],
            [3,2,12],
            [4,3,13],
            [0,4,14],
            [1,10,11],
            [2,11,12],
            [3,12,13],
            [4,13,14],
            [0,14,10],
            [5,6,7],
            [7,5,8],
            [8,5,9],
            [10,11,12],
            [12,10,13],
            [10,13,14]
            ],dtype=int)
        colors=list(idcols)*2+[bkcol]*16   
        scene.visuals.Mesh.__init__(self,vertices=V,faces=I,face_colors=colors,parent=view.scene,shading='flat')
        self.transform=MatrixTransform(transform)


keys=['123QWEASDZXC','456RTYFGHVBN']
axes=np.array([
    [phi,0,1],
    [1,phi,0],
    [0,1,phi],
    [phi,0,-1],
    [-1,phi,0],
    [0,-1,phi],
    [0,1,-phi],
    [-phi,0,1],
    [1,-phi,0],
    [-phi,0,-1],
    [-1,-phi,0],
    [0,-1,-phi]])

clist=np.array([(0,1,1),(0,1,0),(0,0,1),(1,0,1),(1,1,0),(1,0,0)])
class RubiksDodecahedron():
    def __init__(self,turn_speed=.2):
        self.BlockList=[
            TriangleBlock(idcols=clist[tr.tcidlist[i]],transform=tr.ttransformlist[i])
        for i in range(20)]+[ PentagonBlock(idcols=clist[tr.pcidlist[i]],transform=tr.ptransformlist[i]) for i in range(12)]
        self.dirs=np.array([tr.ttransformlist[i].T@np.array([1,1,1,1]) for i in range(20)]+[tr.ptransformlist[i].T@np.array([1,0,phi-1,1]) for i in range(12)])[:,:3]
        self.final_time=None
        self.final_transforms=np.array(tr.ttransformlist+tr.ptransformlist)
        self.axis=None
        self.crotate_objects=None
        self.turn_speed=turn_speed
        self.controls=False
        self.text =[Text(text=keys[0][i]+'/'+keys[1][i],color='white', pos=axes[i]*(phi-1),parent=view.scene,font_size=100,depth_test=True) for i in range(12)]
        for i in range(12):
            self.text[i].visible=False
 

    def begin_turn(self,axis, n_twists):
        self.final_time=time.time()+self.turn_speed*abs(n_twists)
        self.axis=np.sign(n_twists)*axis
        rm=tr.rotation_matrix(axis,n_twists*2*np.pi/5)
        self.crotate_objects=np.nonzero(self.dirs@axis>0)[0]
        self.dirs[self.crotate_objects]=self.dirs[self.crotate_objects]@rm[:3,:3]
        self.final_transforms[self.crotate_objects]=self.final_transforms[self.crotate_objects]@rm

            
    def redraw(self):
        if self.final_time is None:
            pass
        elif time.time()>=self.final_time:
            self.final_time=None
            for i in self.crotate_objects:
                self.BlockList[i].transform=MatrixTransform(self.final_transforms[i])

        else:
            theta=(2*np.pi/5)*(self.final_time-time.time())/self.turn_speed
            rm=tr.rotation_matrix(self.axis,-theta)
            
            for i in self.crotate_objects:
                self.BlockList[i].transform=MatrixTransform(self.final_transforms[i]@rm)
        


    def toggle_controls(self):
        if self.controls==False:
            self.controls=True
            for i in range(12):
                self.text[i].visible=True
        else:
            
            self.controls=False
            for i in range(12):
                self.text[i].visible=False
Dodec=RubiksDodecahedron()

def update(ev):
    Dodec.redraw()    
    canvas.update()



timer = app.Timer()
timer.connect(update)
timer.start(0)



@canvas.events.key_press.connect
def on_key_press(event):
    if (event.key is not None) and (Dodec.final_time is None):
        if event.key==" ":
            Dodec.toggle_controls()
        ind = keys[0].find(str(event.key.name))
        if ind!=-1:
            Dodec.begin_turn(axes[ind],1)
        ind = keys[1].find(str(event.key.name))
        if ind!=-1:
            Dodec.begin_turn(axes[ind],-1)
        

    

if __name__ == '__main__' and sys.flags.interactive == 0:
    app.run()

    
