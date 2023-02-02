import pygame as pg
from numpy import format_float_scientific, sin,cos,pi
import os
import sys
from constants import *



def load_image(image, scale=1):
    fullname = os.path.join("./imgs/", image)
    image = pg.image.load(fullname)

    # size = image.get_size()
    # size = (size[0] * scale, size[1] * scale)
    # image = pg.transform.scale(image, size)
    # image = image.convert()
    
    return image

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

class Pendulum:
    def __init__(self,length,mass,dampcoef,gravity,theta,phi,color="#000000",image = "bitmap.png"):
        self.length = length
        self.w0 = gravity/length
        self.gamma = (length*dampcoef)/mass
        self.theta = theta
        self.phi = phi
        self.color = color
        self.h = 0.0250
        self.image = load_image(image)
    def Auxilaryfun(self,theta,phi):
        return -(self.gamma*phi)-(self.w0*sin(theta))

    def update(self):
        """

        """
        k1 = self.h*self.phi
        l1 = self.h*self.Auxilaryfun(self.theta,self.phi)
        k2 = self.h*(self.phi+(l1*0.5))
        l2 = self.h*self.Auxilaryfun(self.theta+(k1*0.5),self.phi+(l1*0.5))
        k3 = self.h*(self.phi+(l2*0.5))
        l3 = self.h*(self.Auxilaryfun(self.theta+(k2*0.5),self.phi+(l2*0.5)))
        k4 = self.h*(self.phi+l3)
        l4 = self.h*(self.Auxilaryfun(self.theta+k3,self.phi+l3))
        k_ = (1/6)*(k1+k4+2*(k2+k3))
        l_ = (1/6)*(l1+l4+2*(l2+l3))
        self.theta+=k_
        self.phi+=l_

    def draw(self,screen,origin):
        """
        
        """
        x=origin[0]+(self.length*cos((pi*1.5)+self.theta))
        y=origin[1]-(self.length*sin((pi*1.5)+self.theta))
        pg.draw.aaline(screen,self.color,start_pos=origin,end_pos=(x,y))
        screen.blit(self.image,(x-10,y-10))


class PendulumAppro(Pendulum):
    def __init__(self, length, mass, dampcoef, gravity, theta, phi, image,color="#000000"):
        super().__init__(length, mass, dampcoef, gravity, theta, phi,image=image,color=color)

    def Auxilaryfun(self, theta, phi):
        return -((gammal)*phi)-(w0*theta)
    



class DoublePendulum:
    """

    """
    def __init__(self,mass1,mass2,length1,length2,dampcoef,gravity,theta1,theta2,phi1,phi2,color="#000000"):
        self.mass1 = mass1
        self.mass2 = mass2
        self.length1 = length1
        self.length2 = length2
        self.gravity = gravity
        self.theta1 = theta1
        self.theta2 = theta2
        self.phi1 = phi1
        self.phi2 = phi2
        self.h = 0.080
        self.origin = (650,200)
        self.image = load_image("bitmap1.png")
        self.color = color

    def oxillary1(self,theta1,theta2,phi1,phi2):
        diff = theta1-theta2
        return ((-2*sin(diff)*self.mass2*(phi1*phi1*self.length1*cos(diff)+phi2*phi2*self.length2))-self.gravity*(2*self.mass2+self.mass1)*sin(theta1)-self.gravity*self.mass2*sin(theta1-2*theta2))/(self.length1*(2*self.mass1+self.mass2-self.mass2*cos(2*diff)))

    def oxillary2(self,theta1,theta2,phi1,phi2):
        diff = theta1-theta2
        return (2*sin(diff)*(phi1*phi1*self.length1*(self.mass1+self.mass2)+self.gravity*(self.mass1+self.mass2)*cos(theta1)+phi2*phi2*self.length2*self.mass2*cos(diff)))/(self.length2(2*self.mass1+self.mass2-self.mass2*cos(2*diff)))


    def update(self):
        """

        Runge kutta methods
        """
        k11 = self.h*self.phi1
        k12 = self.h*self.phi2
        l11 = self.h*self.oxillary1(self.theta1,self.theta2,self.phi1,self.phi2)
        l12 = self.h*self.oxillary2(self.theta1,self.theta2,self.phi1,self.phi2)
        k21 = self.h*(self.phi1+(l11*0.5))
        k22 = self.h*(self.phi2+(l12*0.5))
        l21 = self.h*self.oxillary1(self.theta1+(k11*0.5),self.theta2+(k12*0.5),self.phi1+(l11*0.5),self.phi2+(l12*0.5))
        l22 = self.h*self.oxillary2(self.theta1+(k11*0.5),self.theta2+(k12*0.5),self.phi1+(l11*0.5),self.phi2+(l12*0.5))
        k31 = self.h*(self.phi1+(l21*0.5))
        k32 = self.h*(self.phi2+(l22*0.5))
        l31 = self.h*self.oxillary1(self.theta1+(k21*0.5),self.theta2+(k22*0.5),self.phi1+(l21*0.5),self.phi2+(l22*0.5))
        l32 = self.h*self.oxillary2(self.theta1+(k21*0.5),self.theta2+(k22*0.5),self.phi1+(l21*0.5),self.phi2+(l22*0.5))
        k41 = self.h*(self.phi1+l31)
        k42 = self.h*(self.phi2+l32)
        l41 = self.h*self.oxillary1(self.theta1+k31,self.theta2+k32,self.phi1+l31,self.phi2+l32)
        l42 = self.h*self.oxillary2(self.theta1+k31,self.theta2+k32,self.phi1+l31,self.phi2+l32)
        k_1 = (1/6)*(k11+k41+2*(k21+k31))
        k_2 = (1/6)*(k21+k42+2*(k22+k32))
        l_1 = (1/6)*(l11+l41+2*(l21+l31))
        l_2 = (1/6)*(l12+l42+2*(l22+l32))

        self.theta1+=k_1
        self.theta2+=k_2
        self.phi1+=l_1
        self.phi2+=l_2
        

    def draw(self,window):
        x1 = self.origin[0]+self.length1*cos((pi*1.5)-self.theta1)
        y1 = self.origin[1]-self.length1*sin((pi*1.5)-self.theta1)
        x2 = x1+self.length2*cos((pi*1.5)-self.theta2)
        y2 =y1-self.length2*sin((pi*1.5)-self.theta2)
        pg.draw.line(window,self.color,start_pos=self.origin,end_pos=(x1+10,y1+10),width =2)
        pg.draw.line(window,self.color,start_pos=(x1+10,y1+10),end_pos=(x2+10,y2+10),width=2)
        window.blit(self.image,(x1,y1))
        window.blit(self.image,(x2,y2))


class Simulation:
    def __init__(self):
        pg.init()
        self.window = pg.display.set_mode((1360,720),pg.RESIZABLE)
        self.size =self.window.get_size()
        self.ff=pg.font.Font("Lato-BoldItalic.ttf",28)
        self.ff2=pg.font.Font("Lato-BoldItalic.ttf",32)
        
        color1 = "#FCE38A"
        color2 = "#C68B59"
        color3 = "#FBC687"
        color4 = "#402218"
        
        self.fg = color4
        self.bg = color1
        self.special = color3
        self.common = color2

        
        self.length1 = 500
        self.length2 = 500
        self.mass1 = 1000
        self.mass2 = 1000
        self.dampcoef = 0.1
        self.gravity = 980
        self.theta1 = 0.550
        self.theta2 = 0.400
        self.phi1 = 0.0
        self.phi2 = 0.0
    def slider(self,x,y,pos):
        cursor = pg.mouse.get_pos()
        clicked_pos = pg.mouse.get_pressed()
        pg.draw.rect(self.window,self.common,(x,y,200,2))
        pg.draw.rect(self.window,self.common,(pos+x-5,y-15,10,30))
        if x+200>cursor[0]>=x and y+25>=cursor[1]>=y-25:
            if clicked_pos[0]==1:
                pg.draw.rect(self.window,self.special,(pos+x-5,y-15,10,30))
                pos = cursor[0]-x

        return pos


    def mainmenu(self):
        run = True
        clock = pg.time.Clock()
        heading = self.ff2.render("Welcome to simulations",True,self.fg,self.special)
        heading_size = heading.get_size()
        heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)

        save = self.ff2.render("run",True,self.fg,self.special)
        save_size = save.get_size()
        save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

        cancel = self.ff2.render("exit",True,self.fg,self.special)
        cancel_size = cancel.get_size()
        cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0]//2,self.size[1]-100,cancel_size[0]+20,heading_size[1]+10),border_radius=5) 

        option1 = self.ff2.render("Single pendulum",True,self.fg,self.special)
        option1_size = option1.get_size()
        option1_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-50,self.size[1]//2-250,200,200),border_radius=15)
        option21 = self.ff2.render("Double pendulum",True,self.fg,self.special)
        option22 = self.ff2.render("(chaotic system)",True,self.fg,self.special)
        option2_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+50,self.size[1]//2-250,200,200),border_radius=15)
        first = 1
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    sys.exit()

                if event.type == pg.MOUSEBUTTONDOWN:
                    if option1_rect.collidepoint(event.pos):
                        self.menuS()
                    if option2_rect.collidepoint(event.pos):
                        first = 0
                        self.menuD()
                    if save_rect.collidepoint(event.pos) and first:
                        self.menuS()
                    elif save_rect.collidepoint(event.pos) and first==0:
                        self.menuD()
        
            clock.tick(60)
            self.window.fill(self.bg)
            self.size = self.window.get_size()
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2-300+5,self.size[1]//2-145,300,200),border_radius=15)
            option1_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-300,self.size[1]//2-150,300,200),border_radius=15)

            pg.draw.rect(self.window,self.common,(self.size[0]//2+55,self.size[1]//2-145,300,200),border_radius=15)
            option2_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+50,self.size[1]//2-150,300,200),border_radius=15)
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2-heading_size[0]//2+5,55,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2+205,self.size[1]-145,save_size[0]+20,save_size[1]+10),border_radius=5)
            save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200,self.size[1]-150,save_size[0]+20,save_size[1]+10),border_radius=5)
            pg.draw.rect(self.window,self.common,(self.size[0]//2-200-cancel_size[0]+5,self.size[1]-145,cancel_size[0]+20,heading_size[1]+10),border_radius=5) 
            cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0],self.size[1]-150,cancel_size[0]+20,heading_size[1]+10),border_radius=5) 


            self.window.blit(heading,(heading_rect.x+10,heading_rect.y+5))
            self.window.blit(option1,(option1_rect.x+35,option1_rect.y+80))
            self.window.blit(option21,(option2_rect.x+35,option2_rect.y+50))
            self.window.blit(option22,(option2_rect.x+40,option2_rect.y+100))
            self.window.blit(save,(save_rect.x+10,save_rect.y+5))
            self.window.blit(cancel,(cancel_rect.x+10,cancel_rect.y+5))


            pg.display.flip()


        pg.quit()
            

    def menuD(self):
        run = True
        clock = pg.time.Clock()
        
        self.size = self.window.get_size()
        heading = self.ff2.render("Menu",True,self.fg,self.special)
        heading_size = heading.get_size()
        heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)

        save = self.ff2.render("run",True,self.fg,self.special)
        save_size = save.get_size()
        save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

        cancel = self.ff2.render("cancel",True,self.fg,self.special)
        cancel_size = cancel.get_size()
        cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0]//2,self.size[1]-100,cancel_size[0]+20,heading_size[1]+10),border_radius=5) 

        length1 = self.length1
        length2 = self.length2
        mass1 = self.mass1
        mass2 = self.mass2
        dampcoef = self.dampcoef
        gravity = self.gravity
        theta1 = self.theta1
        theta2 = self.theta2
        phi1 = self.phi1
        phi2 = self.phi2

        
        option1 = self.ff.render("Length=   "+str(round(length1)),True,self.fg,self.bg)
        option2 = self.ff.render("Length=   "+str(round(length2)),True,self.fg,self.bg)
        option1_size = option1.get_size()

        option3 = self.ff.render("Mass=   "+str(round(mass1)),True,self.fg,self.bg)
        option4 = self.ff.render("Mass=   "+str(round(mass2)),True,self.fg,self.bg)
        option3_size = option3.get_size()


        option5 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef)),True,self.fg,self.bg)
        option5_size = option5.get_size()
        
        option6 = self.ff.render("Gravity=   "+str(round(gravity)),True,self.fg,self.bg)
        option6_size = option6.get_size()

        option7 = self.ff.render("Initial angular displacement 1=   ",True,self.fg,self.bg)
        option8 = self.ff.render("Initial angular displacement 2=   ",True,self.fg,self.bg)
        option7_size = option7.get_size()
        
        option9 = self.ff.render("Initial angular velocity 1=   ",True,self.fg,self.bg)
        option10 = self.ff.render("Initial angular velocity 2=   ",True,self.fg,self.bg)
        option9_size = option9.get_size()
        
        single = self.ff.render("One system",True,self.fg,self.special)
        single_size = single.get_size()
        double = self.ff.render("Two system",True,self.fg,self.special)
        double_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2,550,300,40),border_top_right_radius=5,border_bottom_right_radius=5)


        option1_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0]//2-500,130,10,40))
        option2_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0]//2+50,130,10,40))
        option5_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option5_size[0]//2-500,230,10,40))
        
        input_rect1 = pg.Rect(option1_rect.x+option5_size[0]+10,280, 200, 50)
        input_rect2 = pg.Rect(option2_rect.x+option5_size[0]+10, 280, 200, 50)
        input_rect3 = pg.Rect(option1_rect.x+option5_size[0]+10, 330, 200, 10)
        input_rect4 = pg.Rect(option2_rect.x+option5_size[0]+10,330, 200, 10)



        

        pos1 = length1/10-1
        pos2 = length2/10-1
        pos3 = (mass1-50)/50
        pos4 = (mass2-50)/50
        pos5 = dampcoef
        pos6 = gravity/100

        active1 = False
        active2 = False
        active3 = False
        active4 = False
        user_text1 = str(theta1)
        user_text2 = str(theta2)
        user_text3 = str(phi1)
        user_text4 = str(phi2)


        color_active = self.special
        color_inactive = self.common

        color1 = color_inactive
        color2 = color_inactive
        color3 = color_inactive
        color4 = color_inactive
        
        


        while run:
            for event in pg.event.get():
                if event.type==pg.QUIT:
                    run = False
                    break

                if event.type==pg.KEYDOWN:
                    if active1:
                        if event.key==pg.K_BACKSPACE:
                            user_text1=user_text1[:-1]
                        else:
                            if event.key==pg.K_RETURN or event.key==pg.K_KP_ENTER:
                                active1=False
                                try:
                                    theta1 = float(user_text1)
                                except:
                                    user_text1 = str(theta1)
                            else:
                                user_text1+=event.unicode
                    elif active2:
                        if event.key==pg.K_BACKSPACE:
                            user_text2=user_text2[:-1]
                        else:
                            if event.key==pg.K_RETURN or event.key==pg.K_KP_ENTER:
                                active2=False
                                try:
                                    theta2 = float(user_text2)
                                except:
                                    user_text2 = str(theta2)
                            else:
                                user_text2+=event.unicode

                    elif active3:
                        if event.key==pg.K_BACKSPACE:
                            user_text3=user_text3[:-1]
                        else:
                            if event.key==pg.K_RETURN or event.key==pg.K_KP_ENTER:
                                active3=False
                                try:
                                    phi1 = float(user_text3)
                                except:
                                    user_text3 = str(phi1)
                            else:
                                user_text3+=event.unicode

                    elif active4:
                        if event.key==pg.K_BACKSPACE:
                            user_text4=user_text4[:-1]
                        else:
                            if event.key==pg.K_RETURN or event.key==pg.K_KP_ENTER:
                                active4=False
                                try:
                                    phi2 = float(user_text4)
                                except:
                                    user_text4 = str(phi2)
                            else:
                                user_text4+=event.unicode

                    else:
                        user_text = "0.0"
                if event.type==pg.MOUSEBUTTONDOWN:
                    if input_rect1.collidepoint(event.pos):
                        active1 = True
                        active2 = False
                        active3 = False
                        active4 = False
                    if input_rect2.collidepoint(event.pos):
                        active1 = False
                        active2 = True
                        active3 = False
                        active4 = False
                    if input_rect3.collidepoint(event.pos):
                        active1 = False
                        active2 = False
                        active3 = True
                        active4 = False
                    if input_rect4.collidepoint(event.pos):
                        active1 = False
                        active2 = False
                        active3 = False
                        active4 = True
                    
                    if save_rect.collidepoint(event.pos):
                        self.theta1 = theta1
                        self.theta2 = theta2
                        self.phi1 = phi1
                        self.phi2 = phi2
                        self.runD1()
            clock.tick(60)
            self.window.fill(self.bg)
            self.size = self.window.get_size()

            if active1:
                color1 = color_active
                color2 = color_inactive
                color3 = color_inactive
                color4 = color_inactive
                text_surface1 = self.ff.render(user_text1, True, self.fg,self.special)
                text_surface2 = self.ff.render(user_text2, True, self.fg,self.common)
                text_surface3 = self.ff.render(user_text3, True, self.fg,self.common)
                text_surface4 = self.ff.render(user_text4, True, self.fg,self.common)
            elif active2:
                color1 = color_inactive
                color2 = color_active
                color3 = color_inactive
                color4 = color_inactive
                text_surface1 = self.ff.render(user_text1, True, self.fg,self.common)
                text_surface2 = self.ff.render(user_text2, True, self.fg,self.special)
                text_surface3 = self.ff.render(user_text3, True, self.fg,self.common)
                text_surface4 = self.ff.render(user_text4, True, self.fg,self.common)
            elif active3:
                color1 = color_inactive
                color2 = color_inactive
                color3 = color_active
                color4 = color_inactive
                text_surface1 = self.ff.render(user_text1, True, self.fg,self.common)
                text_surface2 = self.ff.render(user_text2, True, self.fg,self.common)
                text_surface3 = self.ff.render(user_text3, True, self.fg,self.special)
                text_surface4 = self.ff.render(user_text4, True, self.fg,self.common)
            elif active4:
                color1 = color_inactive
                color2 = color_inactive
                color3 = color_inactive
                color4 = color_active
                text_surface1 = self.ff.render(user_text1, True, self.fg,self.common)
                text_surface2 = self.ff.render(user_text2, True, self.fg,self.common)
                text_surface3 = self.ff.render(user_text3, True, self.fg,self.common)
                text_surface4 = self.ff.render(user_text4, True, self.fg,self.special)
            else:
                color1 = color_inactive
                color2 = color_inactive
                color3 = color_inactive
                color4 = color_inactive
                text_surface1 = self.ff.render(user_text1, True, self.fg,self.common)
                text_surface2 = self.ff.render(user_text2, True, self.fg,self.common)
                text_surface3 = self.ff.render(user_text3, True, self.fg,self.common)
                text_surface4 = self.ff.render(user_text4, True, self.fg,self.common)




            input_rect1 = pg.Rect(self.size[0]//2,280, 200, 40)
            input_rect2 = pg.Rect(self.size[0]//2, 330, 200, 40)
            input_rect3 = pg.Rect(self.size[0]//2, 380, 200, 40)
            input_rect4 = pg.Rect(self.size[0]//2,430, 200, 40)
            

            
            
            


            


                
            
            # rectangle drawing
            pg.draw.rect(self.window,self.common,(self.size[0]//2-heading_size[0]//2+5,50+5,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            
            option1_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0]//2-500,130,10,40))
            pos1 = self.slider(option1_rect.x+option5_size[0]+10,130+20,pos1)
            length1 = round((pos1*10)+10)
            option1 = self.ff.render("Length =   "+str(round(length1)),True,self.fg,self.bg)

            option2_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0]//2+50,130,10,40))
            pos2 = self.slider(option2_rect.x+option5_size[0]+10,130+20,pos2)
            length2 = round((pos2*10)+10)
            option2 = self.ff.render("Length =   "+str(round(length2)),True,self.fg,self.bg)
            
            option3_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option3_size[0]//2-500,180,10,40))
            pos3 = self.slider(option1_rect.x+option5_size[0]+10,180+20,pos3)
            mass1 = round(pos3*50+50)
            option3 = self.ff.render("Mass =   "+str(round(mass1)),True,self.fg,self.bg)
            
            option4_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option3_size[0]//2+50,180,10,40))
            pos4 = self.slider(option2_rect.x+option5_size[0]+10,180+20,pos4)
            mass2 = round(pos4*50+50)
            option4 = self.ff.render("Mass =   "+str(round(mass2)),True,self.fg,self.bg)
            
            option5_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option5_size[0]//2-500,230,10,40))
            pos5 = self.slider(option1_rect.x+option5_size[0]+10,230+20,pos5)
            dampcoef = round(pos5/100,3)
            option5 = self.ff.render("Damping Coeffiecient =   "+str(round(dampcoef,3)),True,self.fg,self.bg)
            
            option6_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option6_size[0]//2+50,230,10,40))
            pos6 = self.slider(option2_rect.x+option5_size[0]+10,230+20,pos6)
            gravity = round(100*pos6)
            option6 = self.ff.render("Gravity =   "+str(round(gravity)),True,self.fg,self.bg)

            
            option7_rect = pg.draw.rect(self.window,self.bg,(100,280,10,40))
            option8_rect = pg.draw.rect(self.window,self.bg,(100,330,10,40))
            option9_rect = pg.draw.rect(self.window,self.bg,(100,380,10,40))
            option10_rect = pg.draw.rect(self.window,self.bg,(100,430,10,40))
            


            
            
            



            single_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-single_size[0]-50,500,300,40),border_top_left_radius=5,border_bottom_left_radius=5)
            
            double_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2,500,300,40),border_top_right_radius=5,border_bottom_right_radius=5)
            
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2-200-cancel_size[0]//2+5,self.size[1]-100+5,cancel_size[0]+20,heading_size[1]+10),border_radius=5)
            cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0]//2,self.size[1]-100,cancel_size[0]+20,heading_size[1]+10),border_radius=5)
            pg.draw.rect(self.window,self.common,(self.size[0]//2+200-save_size[0]//2+5,self.size[1]-100+5,save_size[0]+20,save_size[1]+10),border_radius=5)
            save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

       


            
            # blitting
            self.window.blit(heading,(heading_rect.x+10,heading_rect.y+5))
            self.window.blit(option1,option1_rect)
            self.window.blit(option2,option2_rect)
            self.window.blit(option3,option3_rect)
            self.window.blit(option4,option4_rect)
            self.window.blit(option5,option5_rect)
            self.window.blit(option6,option6_rect)
            self.window.blit(option7,option7_rect)
            self.window.blit(option8,option8_rect)
            self.window.blit(option9,option9_rect)
            self.window.blit(option10,option10_rect)
            self.window.blit(single,(single_rect.x+30,single_rect.y+3))
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2,500,5,40)) # line between one pendulum and two pendulum chooser
            self.window.blit(double,(double_rect.x+40,double_rect.y+3))
            self.window.blit(cancel,(cancel_rect.x+10,cancel_rect.y+5))
            self.window.blit(save,(save_rect.x+10,save_rect.y+5))


            input_rect1.w = max(200, text_surface1.get_width()+10)
            input_rect2.w = max(200, text_surface2.get_width()+10)
            input_rect3.w = max(200, text_surface3.get_width()+10)
            input_rect4.w = max(200, text_surface4.get_width()+10)
            pg.draw.rect(self.window, color1, input_rect1)
            pg.draw.rect(self.window, color2, input_rect2)
            pg.draw.rect(self.window, color3, input_rect3)
            pg.draw.rect(self.window, color4, input_rect4)
            self.window.blit(text_surface1, input_rect1)
            self.window.blit(text_surface2, input_rect2)
            self.window.blit(text_surface3, input_rect3)
            self.window.blit(text_surface4, input_rect4)

            pg.display.flip()













        
    def menuS(self):
        run = True
        clock = pg.time.Clock()
        user_text = '0'

        self.size = self.window.get_size()
        heading = self.ff2.render("Menu",True,self.fg,self.special)
        heading_size = heading.get_size()
        heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)

        
        


        save = self.ff2.render("run",True,self.fg,self.special)
        save_size = save.get_size()
        save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

        cancel = self.ff2.render("cancel",True,self.fg,self.special)
        cancel_size = cancel.get_size()
        cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0]//2,self.size[1]-100,cancel_size[0]+20,heading_size[1]+10),border_radius=5) 

        length = self.length1
        mass = self.mass1
        dampcoef = self.dampcoef
        gravity = self.gravity
        theta = self.theta1
        phi = self.phi1

        
        option1 = self.ff.render("Length=   "+str(round(length)),True,self.fg,self.bg)
        option1_size = option1.get_size()

        option2 = self.ff.render("Mass=   "+str(round(mass)),True,self.fg,self.bg)
        option2_size = option2.get_size()


        option3 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef)),True,self.fg,self.bg)
        option3_size = option3.get_size()
        
        option4 = self.ff.render("Gravity=   "+str(round(gravity)),True,self.fg,self.bg)
        option4_size = option4.get_size()

        option5 = self.ff.render("Initial angular displacement=   "+str(theta),True,self.fg,self.bg)
        option5_size = option5.get_size()
        
        option6 = self.ff.render("Initial angular velocity=   "+str(phi),True,self.fg,self.bg)
        option6_size = option6.get_size()
        
        single = self.ff.render("one system",True,self.fg,self.special)
        single_size = single.get_size()
        double = self.ff.render("two system",True,self.fg,self.special)
        double_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2,550,300,40),border_top_right_radius=5,border_bottom_right_radius=5)

        pos1 = length/10-1
        pos2 = (mass-50)/50
        pos3 = dampcoef
        pos4 = gravity/100
        pos5 = ((theta+pi)/pi)*100
        pos6 = phi/10
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                
                if event.type == pg.MOUSEBUTTONDOWN:
                    if save_rect.collidepoint(event.pos):
                        self.length1 = length
                        self.mass1 = mass
                        self.dampcoef = dampcoef
                        self.gravity = gravity
                        self.theta1 = theta
                        self.phi1 = phi
                        self.run1()
                    if double_rect.collidepoint(event.pos):
                        self.menu_of_two_pendulum()
                        
                # all the keys and control. 
            # body
            clock.tick(60)
            self.window.fill(self.bg)
            self.size = self.window.get_size()

            # rectangle drawing
            pg.draw.rect(self.window,self.common,(self.size[0]//2-heading_size[0]//2+5,50+5,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            
            option1_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0]//2+20,130,10,40))
            option2_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option2_size[0]//2+20,180,10,40))
            option3_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option3_size[0]//2-80,230,10,40))
            option4_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option4_size[0]//2-95,280,10,40))
            option5_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option5_size[0]//2,330,10,40))
            option6_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option6_size[0]//2-100,380,10,40))
            single_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-single_size[0]-50,500,300,40),border_top_left_radius=5,border_bottom_left_radius=5)
            
            double_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2,500,300,40),border_top_right_radius=5,border_bottom_right_radius=5)
            
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2-200-cancel_size[0]//2+5,self.size[1]-100+5,cancel_size[0]+20,heading_size[1]+10),border_radius=5)
            cancel_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-cancel_size[0]//2,self.size[1]-100,cancel_size[0]+20,heading_size[1]+10),border_radius=5)
            pg.draw.rect(self.window,self.common,(self.size[0]//2+200-save_size[0]//2+5,self.size[1]-100+5,save_size[0]+20,save_size[1]+10),border_radius=5)
            save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

       


            pos1 = self.slider(self.size[0]//2+option1_size[0],130+20,pos1)
            length = round((pos1*10)+10)

            pos2 = self.slider(self.size[0]//2+option1_size[0],180+20,pos2)
            mass = round(pos2*50+50)

            pos3 = self.slider(self.size[0]//2+option1_size[0],230+20,pos3)
            dampcoef = round(pos3/100,3)
            

            pos4 = self.slider(self.size[0]//2+option1_size[0],280+20,pos4)
            gravity = round(100*pos4)

            pos5 = self.slider(self.size[0]//2+option1_size[0],330+20,pos5)
            theta = round((pos5/100)*pi-pi,5)

            pos6 = self.slider(self.size[0]//2+option1_size[0],380+20,pos6)
            phi = round((pos6/40))
            

            option1 = self.ff.render("Length =   "+str(round(length)),True,self.fg,self.bg)
            option2 = self.ff.render("Mass =   "+str(round(mass)),True,self.fg,self.bg)
            option3 = self.ff.render("Damping Coeffiecient =   "+str(round(dampcoef,3)),True,self.fg,self.bg)
            option4 = self.ff.render("Gravity (in CGS) =   "+str(round(gravity)),True,self.fg,self.bg)
            option5 = self.ff.render("Initial angular displacement =   "+str(round(theta,5)),True,self.fg,self.bg)
            option6 = self.ff.render("Initial angular velocity =   "+str(round(phi)),True,self.fg,self.bg)
            
            # blitting
            self.window.blit(heading,(heading_rect.x+10,heading_rect.y+5))
            self.window.blit(option1,option1_rect)
            self.window.blit(option2,option2_rect)
            self.window.blit(option3,option3_rect)
            self.window.blit(option4,option4_rect)
            self.window.blit(option5,option5_rect)
            self.window.blit(option6,option6_rect)
            self.window.blit(single,(single_rect.x+30,single_rect.y+3))
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2,500,5,40)) # line between one pendulum and two pendulum chooser
            self.window.blit(double,(double_rect.x+40,double_rect.y+3))
            self.window.blit(cancel,(cancel_rect.x+10,cancel_rect.y+5))
            self.window.blit(save,(save_rect.x+10,save_rect.y+5))

            pg.display.flip()
            # pygame quit 


    def menu_of_two_pendulum(self):
        run = True
        clock = pg.time.Clock()
        user_text = '0'

        self.size = self.window.get_size()
        heading = self.ff2.render("Menu for Two Pendulums",True,self.fg,self.special)
        heading_size = heading.get_size()


        save = self.ff2.render("save",True,self.fg,self.special)
        save_size = save.get_size()
        save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+200-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)

        exit = self.ff2.render("X",True,self.fg,self.special)
        exit_size = save.get_size()
        exit_rect = pg.draw.rect(self.window,self.special,(self.size[0],50,save_size[0]+20,save_size[1]+10),border_radius=5)


        length1 = self.length1
        length2 = self.length2
        mass1 = self.mass1
        mass2 = self.mass2
        dampcoef1 = self.dampcoef
        dampcoef2 = self.dampcoef
        
        option11 = self.ff.render("Length=   "+str(round(length1)),True,self.fg,self.bg)
        option12 = self.ff.render("Length=   "+str(round(length2)),True,self.fg,self.bg)
        option1_size = option11.get_size()

        option21 = self.ff.render("Mass=   "+str(round(mass1)),True,self.fg,self.bg)
        option22 = self.ff.render("Mass=   "+str(round(mass2)),True,self.fg,self.bg)
        option2_size = option21.get_size()


        option31 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef1)),True,self.fg,self.bg)
        option32 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef2)),True,self.fg,self.bg)
        option3_size = option31.get_size()


        
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key ==pg.K_ESCAPE:
                        run =False
                        break

                if event.type == pg.MOUSEBUTTONDOWN:
                    if save_rect.collidepoint(event.pos):
                        run = False
                        break
                    
            # here

            clock.tick(60)
            self.window.fill(self.bg)


            heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)
            
            option11_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0],150,80,40))
            option21_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option2_size[0],250,80,40))
            option31_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option3_size[0],350,80,40))

            
            option12_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option1_size[0],150,80,40))
            option22_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option2_size[0],250,80,40))
            option32_rect = pg.draw.rect(self.window,self.bg,(self.size[0]//2-option3_size[0],350,80,40))
            
            save_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-save_size[0]//2,self.size[1]-100,save_size[0]+20,save_size[1]+10),border_radius=5)
            exit_rect = pg.draw.rect(self.window,self.special,(self.size[0],50,save_size[0]+20,save_size[1]+10),border_radius=5)




            option11 = self.ff.render("Length=   "+str(round(length1)),True,self.fg,self.bg)
            option12 = self.ff.render("Length=   "+str(round(length2)),True,self.fg,self.bg)

            option21 = self.ff.render("Mass=   "+str(round(mass1)),True,self.fg,self.bg)
            option22 = self.ff.render("Mass=   "+str(round(mass2)),True,self.fg,self.bg)

            option31 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef1)),True,self.fg,self.bg)
            option32 = self.ff.render("Damping coeffiecient=   "+str(round(dampcoef2)),True,self.fg,self.bg)

            self.window.blit(heading,(heading_rect.x+10,heading_rect.y+5))
            self.window.blit(option11,option11_rect)
            self.window.blit(option12,option12_rect)
            self.window.blit(option21,option21_rect)
            self.window.blit(option22,option22_rect)
            self.window.blit(option31,option31_rect)
            self.window.blit(option32,option32_rect)
            self.window.blit(save,(save_rect.x+10,save_rect.y+5))


            
            pg.display.flip()
        
        
    def run1(self):
        run = True
        clock = pg.time.Clock()
        pen = Pendulum(self.length1,self.mass1,self.dampcoef,self.gravity,self.theta1,self.phi1,image="bitmap1.png")


        menu_option = self.ff2.render("Menu",True,self.fg,self.special)
        exit_option = self.ff2.render("Exit",True,self.fg,self.special)
        length_option = self.ff.render("length = "+str(self.length1)+"cm",True,self.fg,self.bg)
        mass_option = self.ff.render("mass = "+str(self.mass1)+"gm",True,self.fg,self.bg)
        damping_option = self.ff.render("Damping Coeff = "+str(self.dampcoef)+"dyne s/cm",True,self.fg,self.bg)
        
        exit_option = self.ff2.render("Exit",True,self.fg,self.special)
        menu_option_size = menu_option.get_size()
        menu_option_rect = pg.draw.rect(self.window,self.special,(100,100,200,menu_option_size[1]+200))
        exit_option_rect = pg.draw.rect(self.window,self.special,(100,200,200,menu_option_size[1]+20))
        length_option_rect = pg.draw.rect(self.window,self.special,(100,400,menu_option_size[0]+50,menu_option_size[1]+20))
        mass_option_rect = pg.draw.rect(self.window,self.special,(100,500,menu_option_size[0]+50,menu_option_size[1]+20))
        damping_option_rect = pg.draw.rect(self.window,self.special,(100,600,menu_option_size[0]+50,menu_option_size[1]+20))
        
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            pen.draw(self.window,origin)
            pen.update()
            pg.draw.rect(self.window,self.special,(100+5,100+5,200,menu_option_size[1]+20),border_radius=5)
            self.window.blit(menu_option,menu_option_rect)
            self.window.blit(exit_option,exit_option_rect)
            self.window.blit(length_option,length_option_rect)
            self.window.blit(mass_option,mass_option_rect)
            self.window.blit(damping_option,damping_option_rect)
            pg.display.flip()

    def runD1(self):
        run = True
        clock = pg.time.Clock()
        pen = DoublePendulum(self.mass1,self.mass2,self.length1,self.length2,self.dampcoef,self.gravity,self.theta1,self.theta2,self.phi1,self.phi2)
        
        menu_option = self.ff2.render("Menu",True,self.fg,self.special)
        exit_option = self.ff2.render("Exit",True,self.fg,self.special)
        length_option = self.ff.render("length = "+str(self.length1)+"cm",True,self.fg,self.bg)
        mass_option = self.ff.render("mass = "+str(self.mass1)+"gm",True,self.fg,self.bg)
        damping_option = self.ff.render("Damping Coeff = "+str(self.dampcoef)+"dyne s/cm",True,self.fg,self.bg)
        
        exit_option = self.ff2.render("Exit",True,self.fg,self.special)
        menu_option_size = menu_option.get_size()
        menu_option_rect = pg.draw.rect(self.window,self.special,(100,100,200,menu_option_size[1]+200))
        exit_option_rect = pg.draw.rect(self.window,self.special,(100,200,200,menu_option_size[1]+20))
        length_option_rect = pg.draw.rect(self.window,self.special,(100,400,menu_option_size[0]+50,menu_option_size[1]+20))
        mass_option_rect = pg.draw.rect(self.window,self.special,(100,500,menu_option_size[0]+50,menu_option_size[1]+20))
        damping_option_rect = pg.draw.rect(self.window,self.special,(100,600,menu_option_size[0]+50,menu_option_size[1]+20))
        
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            pen.draw(self.window)
            pen.update()
            pg.draw.rect(self.window,self.special,(100+5,100+5,200,menu_option_size[1]+20),border_radius=5)
            self.window.blit(menu_option,menu_option_rect)
            self.window.blit(exit_option,exit_option_rect)
            self.window.blit(length_option,length_option_rect)
            self.window.blit(mass_option,mass_option_rect)
            self.window.blit(damping_option,damping_option_rect)
            pg.display.flip()

    def run2(self,dampcoef1,dampcoef2,gravity1,gravity2):
        run = True
        clock = pg.time.Clock()
        
        pen1 = Pendulum(self.length1,self.mass1,dampcoef1,gravity1,self.theta1,self.phi1,color="#000000",image="bitmap1.png")
        pen2 = PendulumAppro(self.length2,self.mass2,dampcoef2,gravity2,self.theta2,self.phi2,image="bitmap2.png",color="#333333")
        length = "length"
        mass = "mass"
        dampcoeff = "damping"
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
            # main loop 
            clock.tick(60)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            
            pen1.draw(self.window,origin)
            pen2.draw(self.window,origin)
            pen1.update()
            pen2.update()
            pg.display.flip()
        pg.quit()

            
            
                
k = Simulation()
k.mainmenu()
