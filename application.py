import pygame as pg
from numpy import format_float_scientific, sin,cos,pi,sqrt
import os
import sys



def load_image(image, scale=1):
    fullname = os.path.join("./", image)
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
        self.T = ((sqrt(length/gravity))*(1+(0.25*(sin(0.5*theta))**2)+((9/64)*(sin(theta*0.5))**4)))*2*pi
        self.W = (2*pi)/self.T
        self.theta = theta
        self.phi = phi
        self.color = color
        self.mass = mass
        self.gravity = gravity
        self.E0 = 0.5*self.mass*phi*phi+self.mass*self.gravity*self.length*(1-cos(theta))
        self.h = 0.0150
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
    def timeperiod(self):
        return self.T,self.W
    def energy(self):
        T = 0.5*self.mass*self.length*self.length*self.phi*self.phi
        V = self.mass*self.gravity*self.length*(1-cos(self.theta))
        return T,V,T+V
    def initial_E(self):
        return self.E0
        

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
        return -(self.gamma*phi)-(self.w0*theta)


    
    def timeperiod(self):
        W = sqrt(self.gravity/self.length)
        T = (2*pi)/W
        return T,W


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
        self.h = 0.0150
        self.origin = (650,200)
        self.image = load_image("bitmap1.png")
        self.color = color
        
        
        self.l1sq_m1_by2 = 0.5*self.mass1*self.length1**2
        self.l1sq_m2_by2 = 0.5*self.mass2*self.length1**2
        self.l2sq_m2_by2 = 0.5*self.mass2*self.length2**2
        self.l1_l2_m2 = self.length1*self.length2*self.mass2
        T = self.l1sq_m1_by2*self.phi1*self.phi1+self.l1sq_m2_by2*self.phi1*self.phi1+self.l2sq_m2_by2*self.phi2*self.phi2+self.l1_l2_m2*self.phi1*self.phi2*cos(self.theta1-self.theta2)
        V = 2*self.mass1*self.length1*self.gravity+self.mass2*self.length2*self.gravity-(self.mass1+self.mass2)*self.gravity*self.length1*cos(self.theta1)-self.mass2*self.gravity*self.length2*cos(self.theta2)
        self.E0 = T+V


    def oxillary1(self,theta1,theta2,phi1,phi2):
        diff = theta1-theta2
        return (((-2*sin(diff)*self.mass2*(phi1*phi1*self.length1*cos(diff)+phi2*phi2*self.length2))-self.gravity*(2*self.mass2+self.mass1)*sin(theta1)-self.gravity*self.mass2*sin(theta1-2*theta2))/(self.length1*(2*self.mass1+self.mass2-self.mass2*cos(2*diff))))

    def oxillary2(self,theta1,theta2,phi1,phi2):
        diff = theta1-theta2
        return ((2*sin(diff)*(phi1*phi1*self.length1*(self.mass1+self.mass2)+self.gravity*(self.mass1+self.mass2)*cos(theta1)+phi2*phi2*self.length2*self.mass2*cos(diff)))/(self.length2*(2*self.mass1+self.mass2-self.mass2*cos(2*diff))))


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
    def initial_energy(self):
        return self.E0    
    def energy(self):
        T = self.l1sq_m1_by2*self.phi1*self.phi1+self.l1sq_m2_by2*self.phi1*self.phi1+self.l2sq_m2_by2*self.phi2*self.phi2+self.l1_l2_m2*self.phi1*self.phi2*cos(self.theta1-self.theta2)
        V = 2*self.mass1*self.length1*self.gravity+self.mass2*self.length2*self.gravity-(self.mass1+self.mass2)*self.gravity*self.length1*cos(self.theta1)-self.mass2*self.gravity*self.length2*cos(self.theta2)
        E = T+V
        return T,V,E

    def draw(self,window):
        x1 = self.origin[0]+self.length1*cos((pi*1.5)-self.theta1)
        y1 = self.origin[1]-self.length1*sin((pi*1.5)-self.theta1)
        x2 = x1+self.length2*cos((pi*1.5)-self.theta2)
        y2 =y1-self.length2*sin((pi*1.5)-self.theta2)
        pg.draw.aaline(window,self.color,start_pos=self.origin,end_pos=(x1,y1))
        pg.draw.aaline(window,self.color,start_pos=(x1,y1),end_pos=(x2,y2))
        window.blit(self.image,(x1-10,y1-10))
        window.blit(self.image,(x2-10,y2-10))


class Simulation:
    def __init__(self):
        pg.init()
        self.window = pg.display.set_mode((1360,720),pg.RESIZABLE)
        self.size =self.window.get_size()
        pg.display.set_caption("Pendulum simulation")
        self.ff=pg.font.Font("Lato-BoldItalic.ttf",24)
        self.ff2=pg.font.Font("Lato-BoldItalic.ttf",32)
        
        color1 = "#FFF5E4"
        color2 = "#EE6983"
        color3 = "#FFC4C4"
        color4 = "#111111"
        color5 = "#FFC2A2"
        
        self.fg = color4
        self.bg = color1
        self.special = color3
        self.common = color2
        self.extra = color5

        
        self.length1 = 200
        self.lengtho1 = self.length1
        self.length2 = 200
        self.lengtho2 = self.length2
        self.mass1 = 100
        self.masso1 = self.mass1
        self.mass2 = 100
        self.masso2 = self.mass2
        self.dampcoef = 0.0
        self.dampcoefo =self.dampcoef
        self.gravity = 980
        self.gravityo = self.gravity
        self.theta1 = 0.550
        self.thetao1 = self.theta1
        self.theta2 = 0.550
        self.thetao2 = self.theta2
        self.phi1 = 0.0
        self.phio1 = self.phi1
        self.phi2 = 0.0
        self.phio2 = self.phi2
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

    def button_with_shadow(self,text,x,y,background,foreground,shadow_color,font,border_radius=5,shadow_distance=5):
        text = font.render(text,True,foreground,background)
        text_size = text.get_size()
        pg.draw.rect(self.window,shadow_color,(x-text_size[0]//2+shadow_distance,y-text_size[1]//2+shadow_distance,text_size[0]+20,text_size[1]+10),border_radius=border_radius)
        text_rect = pg.draw.rect(self.window,background,(x-text_size[0]//2,y-text_size[1]//2,text_size[0]+20,text_size[1]+10),border_radius=border_radius)
        self.window.blit(text,(text_rect.x+10,text_rect.y+5))
        return text_rect

    def text_left(self,text,x,y,background,foreground,font):
        text = font.render(text,True,foreground,background)
        text_size = text.get_size()
        text_rect = pg.draw.rect(self.window,background,(x,y,text_size[0]+20,text_size[1]+10))
        self.window.blit(text,(text_rect.x+10,text_rect.y+5))

    def back_button(self,size):
        text = self.ff2.render("X",True,self.fg,self.common)
        text_rect = pg.draw.rect(self.window,self.common,(size[0]-100,100,50,50),border_radius=5)
        self.window.blit(text,(text_rect.x+10,text_rect.y+8))
        return text_rect

    def bar(self,x,y,initial_val,new_val):
        pg.draw.rect(self.window,self.special,(x,y,300,15))
        pg.draw.rect(self.window,self.common,(x,y,300*new_val/initial_val,15))

    def inputbox(self,text,input_text,x,y,maximum,activity,background,foreground,active_color,inactive_color,font):
        text = font.render(text,True,foreground,background)
        text_size = text.get_size()
        text_rect = pg.draw.rect(self.window,background,(x,y,text_size[0]+20,text_size[1]+10))
        self.window.blit(text,text_rect)
        if activity:
            text_surface = font.render(input_text, True, foreground,active_color)
            input_rect = pg.Rect(x+text_size[0]+50,y,max(maximum,text_surface.get_width()+10),text_size[1])  
            pg.draw.rect(self.window,active_color,input_rect)
        else:
            text_surface = font.render(input_text, True, foreground,inactive_color)
            input_rect = pg.Rect(x+text_size[0]+50,y,max(maximum,text_surface.get_width()+10),text_size[1])  
            pg.draw.rect(self.window,inactive_color,input_rect)
        self.window.blit(text_surface,input_rect)
        return input_rect
    # def sliderbutton(self,text,x,y,background,foreground,pos,font):
        # text = font.render(text,True,foreground,background)
        # text_size = text.get_size()
        # text_rect = pg.draw.rect(self.window,background,(x,y,text_size[0]+20,text_size[1]+10))
        # self.window.blit(text,(text_rect.x+10,text_rect.y+5))

        # return pos

    
        
    def mainmenu(self):
        run = True
        clock = pg.time.Clock()
        heading_rect = self.button_with_shadow("Welcome to C.M. simulation app",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
        # heading_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-heading_size[0]//2,50,heading_size[0]+20,heading_size[1]+10),border_radius=5)

        exit_rect = self.button_with_shadow("exit",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)

        option1 = self.ff2.render("Single pendulum",True,self.fg,self.special)
        option1_size = option1.get_size()
        option1_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-50,self.size[1]//2-250,200,100),border_radius=15)
        option21 = self.ff2.render("Double pendulum",True,self.fg,self.special)
        option22 = self.ff2.render("(chaotic system)",True,self.fg,self.special)
        option2_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-200-50,self.size[1]//2-250,200,100),border_radius=15)
        option3 = self.ff2.render("Coupled pendulum",True,self.fg,self.special)
        option3_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+50,self.size[1]//2-250,200,100),border_radius=15)
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    sys.exit()

                if event.type == pg.MOUSEBUTTONDOWN:
                    if option1_rect.collidepoint(event.pos):
                        self.menuS()
                    if option2_rect.collidepoint(event.pos):
                        self.menuD()
                    if exit_rect.collidepoint(event.pos):
                        run = False
                        sys.exit()
        
            clock.tick(60)
            self.window.fill(self.bg)
            self.size = self.window.get_size()
            
            pg.draw.rect(self.window,self.common,(self.size[0]//2-300+5,self.size[1]//2-145,300,150),border_radius=15)
            option1_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-300,self.size[1]//2-150,300,150),border_radius=15)

            pg.draw.rect(self.window,self.common,(self.size[0]//2+55,self.size[1]//2-145,300,150),border_radius=15)
            option2_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2+50,self.size[1]//2-150,300,150),border_radius=15)

            pg.draw.rect(self.window,self.common,(self.size[0]//2-300+5,self.size[1]//2+55,300,150),border_radius=15)
            option3_rect = pg.draw.rect(self.window,self.special,(self.size[0]//2-300,self.size[1]//2+50,300,150),border_radius=15)

            heading_rect = self.button_with_shadow("Welcome to C.M. Simulation app",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
            exit_rect = self.button_with_shadow("exit",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)


            self.window.blit(option1,(option1_rect.x+35,option1_rect.y+50))
            self.window.blit(option21,(option2_rect.x+35,option2_rect.y+30))
            self.window.blit(option22,(option2_rect.x+40,option2_rect.y+70))
            self.window.blit(option3,(option3_rect.x+25,option3_rect.y+50))
            pg.display.flip()


        pg.quit()
            
        
    def menuD(self):
        run = True
        clock = pg.time.Clock()
        self.size = self.window.get_size()

        self.length1 = 200
        self.length2 = 200
        self.lengtho1 = 200
        self.lengtho2 = 200
        

        back_button = self.back_button(self.size)
        heading_rect = self.button_with_shadow("Double pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
        back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)


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


        
        pos1 = length1/10-1
        pos2 = length2/10-1
        pos3 = (mass1-20)/20
        pos4 = (mass2-20)/20
        # pos5 = dampcoef*2
        pos5 = gravity/100

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

        
        input_rect1 = self.inputbox("initial displacement 1= ",input_text=str(user_text1),x=100,y=330,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect2 = self.inputbox("initial displacement 2 = ",input_text=str(user_text2),x=100,y=380,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect3 = self.inputbox("initial angular velocity 1 = ",input_text=str(user_text3),x=100,y=430,maximum=200,activity=active3,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect4 = self.inputbox("initial angular velocity 2 = ",input_text=str(user_text4),x=100,y=480,maximum=200,activity=active4,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        one_pen_rect = self.button_with_shadow("one system",self.size[0]//2,self.size[1]-200,self.extra,self.fg,self.common,self.ff,shadow_distance=2)
        two_pen_rect = self.button_with_shadow("two systems",self.size[0]//2+300,self.size[1]-200,self.special,self.fg,self.common,self.ff)

        


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
                    if back_button.collidepoint(event.pos):
                        run = False
                        break
                    if two_pen_rect.collidepoint(event.pos):
                        self.menuD2()
                    if run_rect.collidepoint(event.pos):
                        self.length1 = length1
                        self.length2 = length2
                        self.mass1 = mass1
                        self.mass2 = mass2
                        # self.dampcoef = dampcoef
                        self.gravity = gravity
                        self.theta1 = theta1
                        self.theta2 = theta2
                        self.phi1 = phi1
                        self.phi2 = phi2
                        self.runD1()
                    if back_rect.collidepoint(event.pos):
                        run = False
                        break
                    if reset_rect.collidepoint(event.pos):
                        pos1 = (self.lengtho1/10)-1
                        pos2 = (self.lengtho2/10)-1
                        pos3 = (self.masso1-20)/20
                        pos4 = (self.masso2-20)/20
                        # pos5 = self.dampcoefo*2
                        pos5 = self.gravityo/100
                        user_text1 = str(self.thetao1)
                        user_text2 = str(self.thetao2)
                        user_text3 = str(self.phio1)
                        user_text4 = str(self.phio2)
            clock.tick(120)
            self.window.fill(self.bg)
            self.size = self.window.get_size()

            back_button = self.back_button(self.size)
            heading_rect = self.button_with_shadow("Double pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
            back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            one_pen_rect = self.button_with_shadow("one system",self.size[0]//2+100,self.size[1]-300,self.extra,self.fg,self.common,self.ff,shadow_distance=2)
            two_pen_rect = self.button_with_shadow("two systems",self.size[0]//2+400,self.size[1]-300,self.special,self.fg,self.common,self.ff)


            
            # rectangle drawing
            
            self.text_left("length 1=   "+str(round(length1)),90,170,self.bg,self.fg,self.ff)
            pos1 = self.slider(self.size[0]//2-250,170+20,pos1)
            length1 = round((pos1*10)+10)


            pos2 = self.slider(self.size[0]-250,170+20,pos2)
            length2 = round((pos2*10)+10)
            self.text_left("length 2=   "+str(round(length2)),self.size[0]//2,170,self.bg,self.fg,self.ff)

            pos3 = self.slider(self.size[0]//2-250,220+20,pos3)
            mass1 = round(pos3*20+20)
            self.text_left("mass 1=   "+str(round(mass1)),90,220,self.bg,self.fg,self.ff)

            pos4 = self.slider(self.size[0]-250,220+20,pos4)
            mass2 = round(pos4*20+20)
            self.text_left("mass 2=   "+str(round(mass2)),self.size[0]//2,220,self.bg,self.fg,self.ff)

            # pos5 = self.slider(self.size[0]//2-250,270+20,pos5)
            # dampcoef = round(pos5/200,3)
            # self.text_left("damp coeffciant =   "+str(dampcoef),100,270,self.bg,self.fg,self.ff)

            pos5 = self.slider(self.size[0]//2-250,270+20,pos5)
            gravity = round(10*pos5)+882
            self.text_left("gravity=   "+str(round(gravity)),90,270,self.bg,self.fg,self.ff)

            
            input_rect1 = self.inputbox("initial displacement 1= ",input_text=str(user_text1),x=100,y=330,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect2 = self.inputbox("initial displacement 2 = ",input_text=str(user_text2),x=100,y=380,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect3 = self.inputbox("initial angular velocity 1 = ",input_text=str(user_text3),x=100,y=430,maximum=200,activity=active3,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect4 = self.inputbox("initial angular velocity 2 = ",input_text=str(user_text4),x=100,y=480,maximum=200,activity=active4,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)

            
            pg.display.flip()

    def menuD2(self):
        run = True
        clock = pg.time.Clock()
        self.size = self.window.get_size()
        back_button = self.back_button(self.size)
        while run:
            for event in pg.event.get():
                if event.type==pg.QUIT:
                    run = False
                    break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if back_button.collidepoint(event.pos):
                        run = False
                        break
            
            clock.tick(60)
            self.window.fill(self.bg)
            self.size =self.window.get_size()
            self.text_left("under construction, well don't want to do it, bored ;(",self.size[0]//2-400,self.size[1]//2,self.bg,self.fg,self.ff2)
            back_button = self.back_button(self.size)
            pg.display.flip()
    def menuS(self):
        run = True
        clock = pg.time.Clock()
        self.size = self.window.get_size()
        
        self.length1 = 500
        self.length2 = 500
        self.lengtho1 = 500
        self.lengtho2 = 500
        back_button = self.back_button(self.size)
        heading_rect = self.button_with_shadow("Simple pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
        back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        one_pen_rect = self.button_with_shadow("one pendulum system",self.size[0]//2-200,self.size[1]-200,self.extra,self.fg,self.common,self.ff,shadow_distance=2)
        two_pen_rect = self.button_with_shadow("two pendulum system",self.size[0]//2+200,self.size[1]-200,self.special,self.fg,self.common,self.ff)

        length1 = self.length1
        mass1 = self.mass1
        dampcoef = self.dampcoef
        gravity = self.gravity
        theta1 = self.theta1
        phi1 = self.phi1
        


        pos1 = length1/10-1
        pos2 = (mass1-20)/20
        pos3 = dampcoef*2
        pos4 = gravity/100

        active1 = False
        active2 = False

        user_text1 = str(theta1)
        user_text2 = str(phi1)


        color_active = self.special
        color_inactive = self.common

        
        input_rect1 = self.inputbox("initial displacement= ",input_text=str(user_text1),x=100,y=330,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect2 = self.inputbox("initial angular velocity = ",input_text=str(user_text2),x=100,y=430,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)

            
        
        
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                
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
                                    phi1 = float(user_text2)
                                except:
                                    user_text2 = str(phi1)
                            else:
                                user_text2+=event.unicode
                    else:
                        user_text = "0.0"
                if event.type==pg.MOUSEBUTTONDOWN:
                    if input_rect1.collidepoint(event.pos):
                        active1 = True
                        active2 = False
                    if input_rect2.collidepoint(event.pos):
                        active1 = False
                        active2 = True
                    if back_button.collidepoint(event.pos):
                        run = False
                        break
                    if run_rect.collidepoint(event.pos):
                        self.length1 = length1
                        self.mass1 = mass1
                        self.dampcoef = dampcoef
                        self.gravity = gravity
                        self.theta1 = theta1
                        self.phi1 = phi1
                        self.run1()
                    if back_rect.collidepoint(event.pos):
                        run = False
                        break
                    if reset_rect.collidepoint(event.pos):
                        pos1 = (self.lengtho1/10)-1
                        pos2 = (self.masso1-20)/20
                        pos3 = self.dampcoefo*2
                        pos4 = self.gravityo/100
                        user_text1 = str(self.thetao1)
                        user_text2 = str(self.phio1)
                    if two_pen_rect.collidepoint(event.pos):
                        self.menu_of_two_pendulum()
                        run = False
                        break
                # all the keys and control. 
            # body

            clock.tick(120)
            self.window.fill(self.bg)
            self.size = self.window.get_size()

            back_button = self.back_button(self.size)
            heading_rect = self.button_with_shadow("Simple pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
            back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            one_pen_rect = self.button_with_shadow("one pendulum system",self.size[0]//2-200,self.size[1]-200,self.extra,self.fg,self.common,self.ff,shadow_distance=3)
            two_pen_rect = self.button_with_shadow("two pendulum system",self.size[0]//2+200,self.size[1]-200,self.special,self.fg,self.common,self.ff)

            
            # rectangle drawing
            
            self.text_left("length =   "+str(round(length1)),90,170,self.bg,self.fg,self.ff)
            pos1 = self.slider(self.size[0]//2-250,170+20,pos1)
            length1 = round((pos1*10)+10)


            self.text_left("mass =   "+str(round(mass1)),self.size[0]//2,170,self.bg,self.fg,self.ff)
            pos2 = self.slider(self.size[0]-250,170+20,pos2)
            mass1 = round(pos2*20+20)

            self.text_left("damp coeffciant=   "+str(dampcoef),90,220,self.bg,self.fg,self.ff)
            pos3 = self.slider(self.size[0]//2-250,220+20,pos3)
            dampcoef = round(pos3/200,6)

            self.text_left("gravity=   "+str(round(gravity)),self.size[0]//2,220,self.bg,self.fg,self.ff)
            pos4 = self.slider(self.size[0]-250,220+20,pos4)
            gravity = round(10*pos4)+882

            input_rect1 = self.inputbox("initial displacement= ",input_text=str(user_text1),x=100,y=330,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect2 = self.inputbox("initial angular velocity= ",input_text=str(user_text2),x=100,y=380,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)

            

            pg.display.flip()


    def menu_of_two_pendulum(self):
        run = True
        clock = pg.time.Clock()
        self.size = self.window.get_size()
        
        back_button = self.back_button(self.size)
        heading_rect = self.button_with_shadow("Simple pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
        back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
        exact_rect = self.button_with_shadow("exact solution",self.size[0]//2-200,self.size[1]-200,self.special,self.fg,self.common,self.ff)
        appro_rect = self.button_with_shadow("appoximated solution",self.size[0]//2+200,self.size[1]-200,self.special,self.fg,self.common,self.ff)

        length1 = self.length1
        length2 = self.length2
        mass1 = self.mass1
        mass2 = self.mass2
        dampcoef1 = self.dampcoef
        dampcoef2 = self.dampcoef
        gravity1 = self.gravity
        gravity2 = self.gravity
        theta1 = self.theta1
        theta2 = self.theta2
        phi1 = self.phi1
        phi2 = self.phi2

        
        pos1 = length1/10-1
        pos2 = length2/10-1
        pos3 = (mass1-20)/20
        pos4 = (mass2-20)/20
        pos5 = dampcoef1*2
        pos6 = dampcoef2*2
        pos7 = gravity1/100
        pos8 = gravity2/100


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

        
        input_rect1 = self.inputbox("initial displacement 1= ",input_text=str(user_text1),x=100,y=330,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect2 = self.inputbox("initial displacement 2 = ",input_text=str(user_text2),x=100,y=380,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect3 = self.inputbox("initial angular velocity 1 = ",input_text=str(user_text3),x=100,y=430,maximum=200,activity=active3,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
        input_rect4 = self.inputbox("initial angular velocity 2 = ",input_text=str(user_text4),x=100,y=480,maximum=200,activity=active4,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)

    
        appro = True
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
                    if back_button.collidepoint(event.pos):
                        run = False
                        break
                    if run_rect.collidepoint(event.pos):
                        self.length1 = length1
                        self.length2 = length2
                        self.mass1 = mass1
                        self.mass2 = mass2
                        self.theta1 = theta1
                        self.theta2 = theta2
                        self.phi1 = phi1
                        self.phi2 = phi2
                        if appro:
                            self.run22(dampcoef1,dampcoef2,gravity1,gravity2)
                        else:
                            self.run21(dampcoef1,dampcoef2,gravity1,gravity2)
                    if exact_rect.collidepoint(event.pos):
                        appro=False
                    if appro_rect.collidepoint(event.pos):
                        appro=True
                    if back_rect.collidepoint(event.pos):
                        run = False
                        break
                    if reset_rect.collidepoint(event.pos):
                        pos1 = (self.lengtho1/10)-1
                        pos2 = (self.lengtho2/10)-1
                        pos3 = (self.masso1-20)/20
                        pos4 = (self.masso2-20)/20
                        pos5 = self.dampcoefo*2
                        pos6 = self.dampcoefo*2
                        pos7 = self.gravityo/100
                        pos8 = self.gravityo/100
                        user_text1 = str(self.thetao1)
                        user_text2 = str(self.thetao2)
                        user_text3 = str(self.phio1)
                        user_text4 = str(self.phio2)
            clock.tick(60)
            self.window.fill(self.bg)
            self.size = self.window.get_size()

            back_button = self.back_button(self.size)
            heading_rect = self.button_with_shadow("Simple pendulum",self.size[0]//2,100,self.special,self.fg,self.common,self.ff2)
            back_rect = self.button_with_shadow("back",self.size[0]//2-200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            reset_rect = self.button_with_shadow("reset",self.size[0]//2,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            run_rect = self.button_with_shadow("run",self.size[0]//2+200,self.size[1]-100,self.special,self.fg,self.common,self.ff2)
            if appro:
                exact_rect = self.button_with_shadow("exact solution",self.size[0]//2-200,self.size[1]-180,self.special,self.fg,self.common,self.ff)
                appro_rect = self.button_with_shadow("appoximated solution",self.size[0]//2+200+3,self.size[1]-180+3,self.extra,self.fg,self.common,self.ff,shadow_distance = 3)
            else:
                exact_rect = self.button_with_shadow("exact solution",self.size[0]//2-200+3,self.size[1]-180+3,self.extra,self.fg,self.common,self.ff,shadow_distance=3)
                appro_rect = self.button_with_shadow("appoximated solution",self.size[0]//2+200,self.size[1]-180,self.special,self.fg,self.common,self.ff)
            
            # rectangle drawing
            
            self.text_left("length 1=   "+str(round(length1)),90,170,self.bg,self.fg,self.ff)
            pos1 = self.slider(self.size[0]//2-250,170+20,pos1)
            length1 = round((pos1*10)+10)


            pos2 = self.slider(self.size[0]-250,170+20,pos2)
            length2 = round((pos2*10)+10)
            self.text_left("length 2=   "+str(round(length2)),self.size[0]//2,170,self.bg,self.fg,self.ff)

            pos3 = self.slider(self.size[0]//2-250,220+20,pos3)
            mass1 = round(pos3*20+20)
            self.text_left("mass 1=   "+str(round(mass1)),90,220,self.bg,self.fg,self.ff)

            pos4 = self.slider(self.size[0]-250,220+20,pos4)
            mass2 = round(pos4*20+20)
            self.text_left("mass 2=   "+str(round(mass2)),self.size[0]//2,220,self.bg,self.fg,self.ff)

            pos5 = self.slider(self.size[0]//2-250,270+20,pos5)
            dampcoef1 = round(pos5/200,3)
            self.text_left("damp coeffciant (1)=   "+str(dampcoef1),90,270,self.bg,self.fg,self.ff)

            pos6 = self.slider(self.size[0]-250,270+20,pos6)
            dampcoef2 = round(pos6/200,3)
            self.text_left("damp coeffciant (2)=   "+str(dampcoef2),self.size[0]//2,270,self.bg,self.fg,self.ff)
    
            pos7 = self.slider(self.size[0]//2-250,320+20,pos7)
            gravity1 = round(10*pos7)+882
            self.text_left("gravity (1)=   "+str(round(gravity1)),90,320,self.bg,self.fg,self.ff)

            pos8 = self.slider(self.size[0]-250,320+20,pos8)
            gravity2 = round(10*pos8)+882
            self.text_left("gravity (2)=   "+str(round(gravity2)),self.size[0]//2,320,self.bg,self.fg,self.ff)

            
            input_rect1 = self.inputbox("initial displacement 1= ",input_text=str(user_text1),x=100,y=400,maximum=200,activity=active1,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect2 = self.inputbox("initial displacement 2 = ",input_text=str(user_text2),x=self.size[0]//2,y=400,maximum=200,activity=active2,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect3 = self.inputbox("initial angular velocity 1 = ",input_text=str(user_text3),x=100,y=450,maximum=200,activity=active3,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)
            input_rect4 = self.inputbox("initial angular velocity 2 = ",input_text=str(user_text4),x=self.size[0]//2,y=450,maximum=200,activity=active4,background=self.bg,foreground=self.fg,active_color=color_active,inactive_color=color_inactive,font=self.ff)

            
            pg.display.flip()

        

    def run1(self):
        run = True
        
        clock = pg.time.Clock()
        pen = Pendulum(self.length1,self.mass1,self.dampcoef,self.gravity,self.theta1,self.phi1,image="bitmap1.png")
        E0 = pen.initial_E()
        T,V,E = pen.energy()
        E0 = max(E,E0)
        # T = 2*pi*sqrt(self.length1/self.gravity)
        menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
        TP,W = pen.timeperiod()
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if menu_button.collidepoint(event.pos):
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            
            T,V,E = pen.energy()
            menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
            self.text_left("length= "+str(self.length1)+" cm",50,300,self.bg,self.fg,self.ff)
            self.text_left("mass= "+str(self.mass1)+" gm",50,350,self.bg,self.fg,self.ff)
            self.text_left("time period= "+str(round(TP,4))+" s",50,400,self.bg,self.fg,self.ff)
            self.text_left("angular frequency= "+str(round(W,4))+" rad/s",50,450,self.bg,self.fg,self.ff)
            self.text_left("E= ",50,500,self.bg,self.fg,self.ff)
            self.bar(300,515,E0,E)
            self.text_left("T= ",50,550,self.bg,self.fg,self.ff)
            self.bar(300,565,E0,T)
            self.text_left("V= ",50,600,self.bg,self.fg,self.ff)
            self.bar(300,615,E0,V)
            
            pen.draw(self.window,origin)
            pen.update()

            pg.display.flip()

    



    def run21(self,dampcoef1,dampcoef2,gravity1,gravity2):
        run = True
        clock = pg.time.Clock()
        pen1 = Pendulum(self.length1,self.mass1,dampcoef1,gravity1,self.theta1,self.phi1,image="bitmap1.png")
        pen2 = Pendulum(self.length2,self.mass2,dampcoef2,gravity2,self.theta2,self.phi2,image="bitmap2.png")
        T1,V1,E1 = pen1.energy()
        T2,V2,E2 = pen2.energy()
        E01 = pen1.initial_E()
        E02 = pen2.initial_E()
        E01 = max(E1,E01)
        E02 = max(E2,E02)
        # T = 2*pi*sqrt(self.length1/self.gravity)
        menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
        TP1,W1 = pen1.timeperiod()
        TP2,W2 = pen2.timeperiod()
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if menu_button.collidepoint(event.pos):
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            
            T1,V1,E1 = pen1.energy()
            T2,V2,E2 = pen2.energy()
            menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
            self.text_left("length= "+str(self.length1)+" cm",50,300,self.bg,self.fg,self.ff)
            self.text_left("mass= "+str(self.mass1)+" gm",50,350,self.bg,self.fg,self.ff)
            self.text_left("time period= "+str(round(TP1,4))+" s",50,400,self.bg,self.fg,self.ff)
            self.text_left("angular frequency= "+str(round(W1,4))+" rad/s",50,450,self.bg,self.fg,self.ff)
            self.text_left("E= ",50,500,self.bg,self.fg,self.ff)
            self.bar(100,515,E01,E1)
            self.text_left("T= ",50,550,self.bg,self.fg,self.ff)
            self.bar(100,565,E01,T1)
            self.text_left("V= ",50,600,self.bg,self.fg,self.ff)
            self.bar(100,615,E01,V1)

            self.text_left("length= "+str(self.length2)+" cm",self.size[0]-350,300,self.bg,self.fg,self.ff)
            self.text_left("mass= "+str(self.mass2)+" gm",self.size[0]-350,350,self.bg,self.fg,self.ff)
            self.text_left("time period= "+str(round(TP2,4))+" s",self.size[0]-350,400,self.bg,self.fg,self.ff)
            self.text_left("angular frequency= "+str(round(W2,4))+" rad/s",self.size[0]-350,450,self.bg,self.fg,self.ff)
            # self.text_left("total enetgy= ",self.size[0]-500,500,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-325,515,E02,E2)
            # self.text_left("kinetic energy= ",self.size[0]-500,550,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-325,565,E02,T2)
            # self.text_left("potential energy= ",self.size[0]-500,600,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-325,615,E02,V2)
            
            pen1.draw(self.window,origin)
            pen2.draw(self.window,origin)
            pen1.update()
            pen2.update()

            pg.display.flip()
        
   
    def run22(self,dampcoef1,dampcoef2,gravity1,gravity2):
        run = True
        clock = pg.time.Clock()
        pen1 = Pendulum(self.length1,self.mass1,dampcoef1,gravity1,self.theta1,self.phi1,image="bitmap1.png")
        pen2 = PendulumAppro(self.length2,self.mass2,dampcoef2,gravity2,self.theta2,self.phi2,image="bitmap2.png")
        T1,V1,E1 = pen1.energy()
        T2,V2,E2 = pen2.energy()
        E01 = pen1.initial_E()
        E02 = pen2.initial_E()
        E01 = max(E1,E01)
        E02 = max(E2,E02)
        menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
        TP1,W1 = pen1.timeperiod()
        TP2,W2 = pen2.timeperiod()
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if menu_button.collidepoint(event.pos):
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            
            T1,V1,E1 = pen1.energy()
            T2,V2,E2 = pen2.energy()
            menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
            self.text_left("length= "+str(self.length1)+" cm",50,500,self.bg,self.fg,self.ff)
            self.text_left("mass= "+str(self.mass1)+" gm",50,550,self.bg,self.fg,self.ff)
            self.text_left("time period= "+str(round(TP1,4))+" s",50,600,self.bg,self.fg,self.ff)
            self.text_left("angular frequency= "+str(round(W1,4))+" rad/s",50,650,self.bg,self.fg,self.ff)
            # self.text_left("E= ",50,500,self.bg,self.fg,self.ff)
            # self.bar(100,515,E01,E1)
            # self.text_left("T= ",50,550,self.bg,self.fg,self.ff)
            # self.bar(100,565,E01,T1)
            # self.text_left("V= ",50,600,self.bg,self.fg,self.ff)
            # self.bar(100,615,E01,V1)

            self.text_left("length= "+str(self.length2)+" cm",self.size[0]-350,500,self.bg,self.fg,self.ff)
            self.text_left("mass= "+str(self.mass2)+" gm",self.size[0]-350,550,self.bg,self.fg,self.ff)
            self.text_left("time period= "+str(round(TP2,4))+" s",self.size[0]-350,600,self.bg,self.fg,self.ff)
            self.text_left("angular frequency= "+str(round(W2,4))+" rad/s",self.size[0]-350,650,self.bg,self.fg,self.ff)
            # self.text_left("total enetgy= ",self.size[0]-500,500,self.bg,self.fg,self.ff)
            # self.bar(self.size[0]-315,515,E02,E2)
            # self.text_left("kinetic energy= ",self.size[0]-500,550,self.bg,self.fg,self.ff)
            # self.bar(self.size[0]-315,565,E02,T2)
            # self.text_left("potential energy= ",self.size[0]-500,600,self.bg,self.fg,self.ff)
            # self.bar(self.size[0]-315,615,E02,V2)
            
            pen1.draw(self.window,origin)
            pen2.draw(self.window,origin)
            pen1.update()
            pen2.update()

            pg.display.flip()
        

    def runD1(self):
        run = True
        clock = pg.time.Clock()
        pen = DoublePendulum(self.mass1,self.mass2,self.length1,self.length2,self.dampcoef,self.gravity,self.theta1,self.theta2,self.phi1,self.phi2)
        
        menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
        T,V,E = pen.energy()
        E0 = pen.initial_energy()
        E0 = max(E,E0)
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if menu_button.collidepoint(event.pos):
                        run = False
                        break
            # main loop 
            clock.tick(120)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            
            menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
            T,V,E = pen.energy()
            self.text_left("length1= "+str(self.length1)+" cm",50,400,self.bg,self.fg,self.ff)
            self.text_left("length2= "+str(self.length2)+" cm",50,450,self.bg,self.fg,self.ff)
            self.text_left("mass1= "+str(self.mass1)+" gm",50,500,self.bg,self.fg,self.ff)
            self.text_left("mass2= "+str(self.mass2)+" gm",50,550,self.bg,self.fg,self.ff)
            self.text_left("E",self.size[0]-50,500,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-400,515,E0,E)
            self.text_left("T ",self.size[0]-50,550,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-400,565,E0,T)
            self.text_left("V ",self.size[0]-50,600,self.bg,self.fg,self.ff)
            self.bar(self.size[0]-400,615,E0,V)
            pen.draw(self.window)
            pen.update()

            pg.display.flip()

    def runD2(self,theta11,theta21,phi11,phi21,theta12,theta22,phi12,phi22):
        run = True
        clock = pg.time.Clock()
        pen1 = DoublePendulum(self.mass1,self.mass2,self.length1,self.length2,self.dampcoef,self.gravity,theta11,theta21,phi11,phi21)
        pen2 = DoublePendulum(self.mass1,self.mass2,self.length1,self.length2,self.dampcoef,self.gravity,theta12,theta22,phi12,phi22) 
        
        menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
        # T,V,E = pen1.energy()
        # E0 = pen.initial_energy()
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break
                if event.type == pg.KEYDOWN:
                    if event.key==pg.K_ESCAPE:
                        run = False
                        break
                if event.type==pg.MOUSEBUTTONDOWN:
                    if menu_button.collidepoint(event.pos):
                        run = False
                        break
            # main loop 
            clock.tick(60)
            self.window.fill(self.bg)
            origin = (self.window.get_size()[0]//2,100)
            menu_button = self.button_with_shadow("menu",100,100,self.special,self.fg,self.common,self.ff2)
            # T,V,E = pen.energy()
            self.text_left("length1= "+str(self.length1)+" cm",50,400,self.bg,self.fg,self.ff)
            self.text_left("length2= "+str(self.length2)+" cm",50,450,self.bg,self.fg,self.ff)
            self.text_left("mass1= "+str(self.mass1)+" gm",50,500,self.bg,self.fg,self.ff)
            self.text_left("mass2= "+str(self.mass2)+" gm",50,550,self.bg,self.fg,self.ff)
            # self.text_left("E= ",50,500,self.bg,self.fg,self.ff)
            # self.bar(100,515,E0,E)
            # self.text_left("T= ",50,550,self.bg,self.fg,self.ff)
            # self.bar(100,565,E0,T)
            # self.text_left("V= ",50,600,self.bg,self.fg,self.ff)
            # self.bar(100,615,E0,V)
            pen1.draw(self.window)
            pen2.draw(self.window)
            pen1.update()
            pen2.update()

            pg.display.flip()
    

            
            
                
k = Simulation()
k.mainmenu()
