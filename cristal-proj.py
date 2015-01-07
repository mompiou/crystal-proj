#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 18:03:26 2012

@author: ania-fred
"""


from __future__ import division
import numpy as np
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from Tkinter import *
from tkFileDialog import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import ttk
from mpl_toolkits.mplot3d import Axes3D
import os

def Rot(th,a,b,c):
   
   aa=a/np.linalg.norm([a,b,c]);
   bb=b/np.linalg.norm([a,b,c]);
   cc=c/np.linalg.norm([a,b,c]);
   c1=np.array([[1,0,0],[0,1,0],[0,0,1]],float)
   c2=np.array([[aa**2,aa*bb,aa*cc],[bb*aa,bb**2,bb*cc],[cc*aa,
                cc*bb,cc**2]],float)
   c3=np.array([[0,-cc,bb],[cc,0,-aa],[-bb,aa,0]],float)
   R=np.cos(th)*c1+(1-np.cos(th))*c2+np.sin(th)*c3

   return R    




########################################################"
#(xyz)=tD^-1(hkl)=Dstar(hkl) (1er transfo du plan hkl vers les coordonnees cartesiennes)
#(uvw)=D^-1(xyz)=tDstar(xyz)=(xyz)Dstar (2eme transfo des coord cartesiennes vers les directions)
###################################################  

#coordonnees dans le plan
def coord_ortho(P,plan):
    global Dstar, varname
    ref=np.array([0,0,1])
    planN=np.dot(Dstar,plan)
    the=np.arccos(np.dot(planN,ref)/(np.linalg.norm(ref)*np.linalg.norm(planN)))      
    if plan[0]==0 and plan[1]==0:
        axe=ref
    else:
        axe=np.cross(planN,ref)
        
    M=np.dot(Rot(the,axe[0],axe[1],axe[2]),P)
    #print(the*180/np.pi)
    return M

def unique(a):
    a = np.sort(a)
    b = np.diff(a)
    b = np.r_[1, b]
    return a[b != 0]



def calcul():
  
    global vec,varname,atom0,Dstar,taille,zoom, EL, Dz

   
   
    if varname!=0:
        f_space=open(varname,"r")
    else:
        varname=getFileName()
        f_space=open(varname,"r")
        
    crist=[]
    
    for line in f_space:
        crist.append(map(str, line.split()))
        
           
    f_space.close()    
    vec=[]
    atom0=[]
    for i in range(0, np.size(crist)):
        if np.size(crist[i])==3:
            vec.append(crist[i])
        else:
            atom0.append(crist[i])
            
    vec=np.array(vec,float)
    atom0=np.array(atom0,float)
    
    maxi=np.int(atom0[np.shape(atom0)[0]-1,0])
    E=np.array([0,0,0,0])
    EL=np.array([0,0,0,0,0])
    H=np.array([[0]])
    Dz=np.array([0,0,0])
    
    for h in range(1,maxi+1):
        Ca=calcul_atom(atom0[atom0[:,0] == h])
       
        E=np.vstack((E,Ca[0]))
        Dz=np.vstack((Dz,Ca[1]))        
        H=np.vstack((H,h*np.ones((np.shape(Ca[0])[0],1))))
    
    EL=np.append(E,H,axis=1)
    EL=np.delete(EL,(0),axis=0)  
    Dz=np.delete(Dz,(0),axis=0)  
 


def trace():
    global vec,varname,atom0,Dstar,taille,zoom, EL, Dz
       
    fi = f.add_subplot(111) 
    fi.figure.clear()
    fi = f.add_subplot(111) 
    sim=taille.get()
    
    if rond.get()==1:    
        fi.scatter(EL[:,0],EL[:,1],s=sim,c=EL[:,3],marker='o')
    else:
        fi.scatter(EL[:,0],EL[:,1],s=sim,c=EL[:,3],marker='s')
    
    
   
    if ato.get()==1:    
        for k in range(0,np.shape(EL)[0]):
            fi.annotate(str(int(EL[k,4])),(EL[k,0],EL[k,1]))    
    
    
    if lab.get()==1:
        for q in range(0,np.shape(EL)[0]):
            at=Dz[q,:]
            at=np.dot(at,Dstar)
            
            vector=str(np.around(at[0],decimals=3))+','+str(np.around(at[1],decimals=3))+','+str(np.around(at[2],decimals=3))
            fi.annotate(vector,(EL[q,0],EL[q,1]))
            

    
    fi.axis('off') 
    fi.axis('equal')
    fi.figure.canvas.draw() 

   

     
        
    


def rep():
    global varname, vec,E,C,Dz,atom0
    fi = f.add_subplot(111) 
    fi.figure.clear()
    fi = f.add_subplot(111,projection='3d') 
    sim=taille.get()
    if varname!=0:
        f_space=open(varname,"r")
    else:
        varname=getFileName()
        f_space=open(varname,"r")
        
    crist=[]
    
    for line in f_space:
        crist.append(map(str, line.split()))
        
           
    f_space.close()    
    vec=[]
    atom0=[]
    for i in range(0, np.size(crist)):
        if np.size(crist[i])==3:
            vec.append(crist[i])
        else:
            atom0.append(crist[i])
            
    vec=np.array(vec,float)
    atom0=np.array(atom0,float)
    
    maxi=np.int(atom0[np.shape(atom0)[0]-1,0])
    
    for h in range(1,maxi+1):

        E=calcul_rep(atom0[atom0[:,0] == h])
                
        fi.scatter(E[:,0],E[:,1],E[:,2],s=sim,c=str(h/maxi))
           
      
    
        
    fi.axis('off') 
    fi.figure.canvas.draw() 
    
def calcul_rep(atom):  
    global Dstar,varname,C, D0,Dz,planN,plan,vec,c

    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
       
    alp=alp*np.pi/180;
    bet=bet*np.pi/180;
    gam=gam*np.pi/180;
   
    V=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
    D=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,V/(a*b*np.sin(gam))]])
    Dstar=np.transpose(np.linalg.inv(D))
    
        
    na_rep=eval(size_entry_a.get())
    nb_rep=eval(size_entry_b.get())
    nc_rep=eval(size_entry_c.get())
    
    A=np.zeros((np.shape(atom)[0],np.shape(atom)[1]-1))
    w=0
    
    for v in range(0,np.shape(atom)[0]):
        
        A[w,:]=np.dot(D,np.array([atom[v,1],atom[v,2],atom[v,3]]))
        w=w+1
    
    atom_pos=np.array(A[0,:])
    for f in range(0,np.shape(A)[0]):
        for i in range(-na_rep,na_rep+1):
            for j in range(-nb_rep,nb_rep+1):
                for k in range(-nc_rep,nc_rep+1):
                        
                        atom_pos=np.vstack((atom_pos,A[f,:]+i*a*vec[0,:]+j*b*vec[1,:]+k*c*vec[2,:]))

     
    return atom_pos    
    

    
def calcul_atom(atom):
    global Dstar,varname,C, D0,planN,plan,vec, atom_pos

    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*np.pi/180
    bet=bet*np.pi/180
    gam=gam*np.pi/180
  
   
    V=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
    D=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,V/(a*b*np.sin(gam))]])
    Dstar=np.transpose(np.linalg.inv(D))
    
    #atom positions
    
    
    
    
    na=eval(size_entry_a.get())
    nb=eval(size_entry_b.get())
    nc=eval(size_entry_c.get())
    
    A=np.zeros((np.shape(atom)[0],np.shape(atom)[1]-1))
    w=0
    
    for v in range(0,np.shape(atom)[0]):
        
        A[w,:]=np.dot(D,np.array([atom[v,1],atom[v,2],atom[v,3]]))
        w=w+1
    
    atom_pos=np.array(A[0,:])
#    print(atom_pos)
    for f in range(0,np.shape(A)[0]):
        for i in range(-na,na+1):
            for j in range(-nb,nb+1):
                for k in range(-nc,nc+1):
                        
                        atom_pos=np.vstack((atom_pos,A[f,:]+i*a*vec[0,:]+j*b*vec[1,:]+k*c*vec[2,:]))
    h=eval(h_entry.get())
    k=eval(k_entry.get())
    l=eval(l_entry.get())
    

    plan=np.array([h,k,l])
    planN=np.dot(Dstar,plan)
    #print(planN)
    Dz=np.array([0,0,0])
    D0=np.array([0])
    C=np.array([0,0,0])
    
    L=np.zeros(np.shape(atom_pos)[0])
    tt=0
    for t in range(0,np.shape(atom_pos)[0]):
        L[tt]=np.around(np.dot(planN,atom_pos[t]),decimals=4)
        tt=tt+1
    
    Le=unique(np.abs(L))
    #print(Le)
    cc=eval(n1_entry.get())
    dd=eval(n2_entry.get())
    
    for y in range(cc,dd):
        for i in range(0,np.shape(atom_pos)[0]):
       
           if (np.around(np.dot(planN,atom_pos[i]),decimals=4))==Le[y]:
       
              Dz=np.vstack((Dz,atom_pos[i]))
         
    
        
    for j in range(1,np.shape(Dz)[0]):
        C=np.vstack((C,coord_ortho(Dz[j,:],plan)))
        D0=np.vstack((D0,1+np.abs(np.around(np.dot(planN,Dz[j]),decimals=4))))
    C=np.delete(C,(0),axis=0)    
    D0=np.delete(D0,(0),axis=0)
   
    F=np.append(C,D0,axis=1)
           
    return F, Dz, atom_pos    


def getFileName():
    global varname
    varname = askopenfilename()
    return varname



   
root = Tk()
root.title('Plane')
root.geometry('1102x839')


##################"
f = Figure(facecolor='white',figsize=[2,2],dpi=100)

canvas = FigureCanvasTkAgg(f, master=root)

canvas.get_tk_widget().place(relx=0.33,rely=0.0,relheight=1.0,relwidth=0.67)

canvas.show()
toolbar = NavigationToolbar2TkAgg( canvas, root )
toolbar.zoom('on')
toolbar.update()
######################################


style = ttk.Style()
theme = style.theme_use()
default = style.lookup(theme, 'background')
root.configure(background=default)

def init():
    global varname, lab, ato, rond,zoom,taille
    varname=0
    lab=IntVar()
    ato=IntVar()
    rond=IntVar()
    zoom=IntVar()
    taille=IntVar()
    
init()  

plot_button = Button (master=root)
plot_button.place(relx=0.03,rely=0.35,height=27,width=114)
plot_button.configure(activebackground="#f9f9f9")
plot_button.configure(background="#00ff00")
plot_button.configure(text='''Draw plane''',command=trace)

a_label = Label (master=root)
a_label.place(relx=0.03,rely=0.05,height=19,width=11)
a_label.configure(activebackground="#f9f9f9")
a_label.configure(text='''a''')

b_label = Label (master=root)
b_label.place(relx=0.03,rely=0.09,height=19,width=12)
b_label.configure(activebackground="#f9f9f9")
b_label.configure(text='''b''')

c_label = Label (master=root)
c_label.place(relx=0.03,rely=0.12,height=19,width=11)
c_label.configure(activebackground="#f9f9f9")
c_label.configure(text='''c''')

a_entry = Entry (master=root)
a_entry.place(relx=0.05,rely=0.05,relheight=0.02,relwidth=0.04)
a_entry.configure(selectbackground="#c4c4c4")

b_entry = Entry (master=root)
b_entry.place(relx=0.05,rely=0.09,relheight=0.02,relwidth=0.04)
b_entry.configure(selectbackground="#c4c4c4")

c_entry = Entry (master=root)
c_entry.place(relx=0.05,rely=0.12,relheight=0.02,relwidth=0.04)
c_entry.configure(selectbackground="#c4c4c4")

alp_label = Label (master=root)
alp_label.place(relx=0.02,rely=0.16,height=19,width=37)
alp_label.configure(activebackground="#f9f9f9")
alp_label.configure(text='''alpha''')

beta_label = Label (master=root)
beta_label.place(relx=0.02,rely=0.2,height=19,width=31)
beta_label.configure(activebackground="#f9f9f9")
beta_label.configure(text='''beta''')

gamma_label = Label (master=root)
gamma_label.place(relx=0.01,rely=0.23,height=19,width=52)
gamma_label.configure(activebackground="#f9f9f9")
gamma_label.configure(text='''gamma''')

alp_entry = Entry (master=root)
alp_entry.place(relx=0.05,rely=0.16,relheight=0.02,relwidth=0.04)
alp_entry.configure(selectbackground="#c4c4c4")

bet_entry = Entry (master=root)
bet_entry.place(relx=0.05,rely=0.2,relheight=0.02,relwidth=0.04)
bet_entry.configure(selectbackground="#c4c4c4")

gam_entry = Entry (master=root)
gam_entry.place(relx=0.05,rely=0.23,relheight=0.02,relwidth=0.04)
gam_entry.configure(selectbackground="#c4c4c4")

plan_label = Label (master=root)
plan_label.place(relx=0.15,rely=0.05,height=19,width=61)
plan_label.configure(activebackground="#f9f9f9")
plan_label.configure(text='''Plane (hkl)''')

h_entry = Entry (master=root)
h_entry.place(relx=0.2,rely=0.05,relheight=0.02,relwidth=0.02)
h_entry.configure(selectbackground="#c4c4c4")

k_entry = Entry (master=root)
k_entry.place(relx=0.24,rely=0.05,relheight=0.02,relwidth=0.02)
k_entry.configure(selectbackground="#c4c4c4")

l_entry = Entry (master=root)
l_entry.place(relx=0.28,rely=0.05,relheight=0.02,relwidth=0.02)
l_entry.configure(selectbackground="#c4c4c4")

n_label = Label (master=root)
n_label.place(relx=0.15,rely=0.11,height=19,width=55)
n_label.configure(activebackground="#f9f9f9")
n_label.configure(text='''Layers''')

n1_entry = Entry (master=root)
n1_entry.place(relx=0.21,rely=0.11,relheight=0.02,relwidth=0.02)
n1_entry.configure(selectbackground="#c4c4c4")

n2_entry = Entry (master=root)
n2_entry.place(relx=0.24,rely=0.11,relheight=0.02,relwidth=0.02)
n2_entry.configure(selectbackground="#c4c4c4")

size_label = Label (master=root)
size_label.place(relx=0.13,rely=0.16,height=19,width=94)
size_label.configure(activebackground="#f9f9f9")
size_label.configure(text='''Size''')

size_entry_a = Entry (master=root)
size_entry_a.place(relx=0.21,rely=0.16,relheight=0.02
        ,relwidth=0.02)
size_entry_a.configure(selectbackground="#c4c4c4")


size_entry_b = Entry (master=root)
size_entry_b.place(relx=0.25,rely=0.16,relheight=0.02
        ,relwidth=0.02)
size_entry_b.configure(selectbackground="#c4c4c4")

size_entry_c = Entry (master=root)
size_entry_c.place(relx=0.28,rely=0.16,relheight=0.02
        ,relwidth=0.02)
size_entry_c.configure(selectbackground="#c4c4c4")

calcul_button = Button (master=root)
calcul_button.place(relx=0.03,rely=0.31,height=27,width=114)
calcul_button.configure(activebackground="#ff0000")
calcul_button.configure(background="#ff0000")
calcul_button.configure(text='''Calculate''',command=calcul)

iconButton = Button (master=root)
iconButton.place(relx=0.15,rely=0.23,height=27,width=172)
iconButton.configure(activebackground="#f9f9f9")
iconButton.configure(text='''Change the structure''',command=getFileName)

label_check = Checkbutton (master=root)
label_check.place(relx=0.22,rely=0.36,relheight=0.02,relwidth=0.07)
label_check.configure(text='''Labels''')
label_check.configure(variable=lab)

atom_check = Checkbutton (master=root)
atom_check.place(relx=0.22,rely=0.40,relheight=0.02,relwidth=0.07)
atom_check.configure(activebackground="#f9f9f9")
atom_check.configure(text='''Atoms''')
atom_check.configure(variable=ato)


rep_button = Button (master=root)
rep_button.place(relx=0.03,rely=0.4,height=27,width=134)
rep_button.configure(activebackground="#f9f9f9")
rep_button.configure(background="#00ffff")
rep_button.configure(text='''3D''',command=rep)


#zoom_scale = Scale (master=root,from_=0, to=200, orient=HORIZONTAL,variable=zoom)
#zoom_scale.place(relx=0.10,rely=0.57,relheight=0.06,relwidth=0.07)
#
#
#zoom_label = Label (master=root)
#zoom_label.place(relx=0.04,rely=0.59,height=19,width=38)
#zoom_label.configure(text='''Zoom''')

rond_check = Checkbutton (master=root)
rond_check.place(relx=0.02,rely=0.55,relheight=0.02,relwidth=0.15)
rond_check.configure(activebackground="#f9f9f9")
rond_check.configure(text='''Circles/Squares''')
rond_check.configure(variable=rond)

size_scale = Scale (master=root,from_=0, to=200, orient=HORIZONTAL,variable=taille)
size_scale.place(relx=0.10,rely=0.49,relheight=0.06,relwidth=0.07)
size_scale.configure(activebackground="#f9f9f9")
size_scale.configure(troughcolor="#c4c4c4")

size_marker_label = Label (master=root)
size_marker_label.place(relx=0.01,rely=0.51,height=19,width=104)
size_marker_label.configure(activebackground="#f9f9f9")
size_marker_label.configure(text='''Marker size''')

menu = Menu(master=root)
root.config(menu=menu)

def createstructure(i):
    return lambda:structure(i)    
    
    
def structure(i0):
    global x0
    
    a_entry.delete(0,END)
    a_entry.insert(1,eval(x0[i0][1]))
    b_entry.delete(0,END)    
    b_entry.insert(1,eval(x0[i0][2]))
    c_entry.delete(0,END)    
    c_entry.insert(1,eval(x0[i0][3]))
    alp_entry.delete(0,END)    
    alp_entry.insert(1,eval(x0[i0][4]))
    bet_entry.delete(0,END)    
    bet_entry.insert(1,eval(x0[i0][5]))
    gam_entry.delete(0,END)    
    gam_entry.insert(1,eval(x0[i0][6]))
    
cristalmenu=Menu(menu,tearoff=0)
menu.add_cascade(label="Structures", menu=cristalmenu)
file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")

x0=[]
i=0
for line in file_struct:
    x0.append(map(str, line.split()))
    cristalmenu.add_command(label=x0[i][0], command=createstructure(i))
    i=i+1


##############################################################
# fonction pour quitter
#############################################################
def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#############################################################


 
a_entry.insert(1,1)
b_entry.insert(1,1)
c_entry.insert(1,1)
h_entry.insert(1,1)
k_entry.insert(0,0)
l_entry.insert(0,0)
alp_entry.insert(90,90)
bet_entry.insert(90,90)
gam_entry.insert(90,90)
n1_entry.insert(0,0)
n2_entry.insert(2,2)
size_entry_a.insert(2,2)
size_entry_b.insert(2,2)
size_entry_c.insert(2,2)
size_scale.set(50)
rond_check.select()







root.mainloop()
