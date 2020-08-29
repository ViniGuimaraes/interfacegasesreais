from tkinter import *
from math import *
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

global aux_mudar_constante
aux_mudar_constante = ["Pa", "m^3"]

def GerarIsotermas(tfn_isot, lbrR_isot, tfT_isot, tfa_isot, tfb_isot):
    frToolBarBottom.grid_forget()
    Rv = lbrR_isot.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])

    n = float(tfn_isot.get())
    T = float(tfT_isot.get())
    a = float(tfa_isot.get())
    b = float(tfb_isot.get())

    Pc_vdw = a/(27*b**2)
    Vc_vdw = 3*n*b
    Tc_vdw = (8*a)/(27*R*b)
    Ttio_vdw = T/Tc_vdw

    Vc_b = 3*n*b
    Tc_b = ((8*a)/(27*R*b))**(1/2)
    Pc_b = ((n*R*Tc_b)/(Vc_b-n*b))-((a*n**2)/(Tc_b*Vc_b**2))
    Ttio_b = T/Tc_b

    Pc_d = a/(4*exp(2)*b**2)
    Vc_d = 2*n*b
    Tc_d = a/(4*R*b)
    Ttio_d = T/Tc_d

    Vc_rk = (1+2**(1/3)+2**(2/3))*b*n
    Tc_rk = ((3*a*b*n**2)/(R*Vc_rk**2))**(2/3)
    Pc_rk = (R*Tc_rk)/(3*Vc_rk)
    Ttio_rk = T/Tc_rk

    V = np.linspace(0.5, 4, 100)

    Vtio = np.linspace(0.5, 4, 100)
    
    s0 = ((8*Ttio_vdw)/(3*Vtio-1))-3/Vtio**2
    s1 = (1/Pc_b)*(((n*R*(Tc_b*Ttio_b))/((Vc_b*Vtio)-n*b))-((a*n**2)/((Tc_b*Ttio_b)*(Vc_b*Vtio)**2)))
    s2 = (1/Pc_rk)*(((n*R*(Ttio_rk*Tc_rk))/((Vc_rk*Vtio)-n*b))-((a*n**2)/(((Vc_rk*Vtio)+n*b)*(Vc_rk*Vtio)*T**(1/2))))
    s3 = (1/Pc_d)*(((n*R*(Ttio_d*Tc_d))/((Vtio*Vc_d)-n*b))*(np.exp((-a*n**2)/((Ttio_d*Tc_d)*R*(Vtio*Vc_d)))))

    fig, ax = plt.subplots()
    l0, = ax.plot(Vtio, s0, lw=2, color='k', label='VdW')
    l1, = ax.plot(Vtio, s1, visible=False, lw=2, color='r', label='Berth')
    l2, = ax.plot(Vtio, s2, visible=False, lw=2, color='g', label='RK')
    l3, = ax.plot(Vtio, s3, visible=False, lw=2, color='y', label='Diet')
    plt.subplots_adjust(left=0.3)

    plt.ylim([0, 3])

    plt.title("Isotermas")
    plt.xlabel('$V/V_c$')
    plt.ylabel('$P/P_c$')
    plt.grid(alpha = .4, linestyle = "--")

    lines = [l0, l1, l2, l3]

    # Make checkbuttons with all plotted lines with correct visibility
    rax = plt.axes([0.05, 0.4, 0.1, 0.15])
    labels = [str(line.get_label()) for line in lines]
    visibility = [line.get_visible() for line in lines]
    check = CheckButtons(rax, labels, visibility)
    for i, c in enumerate(["k", "r", "g", "y"]):
        check.labels[i].set_color(c)
        check.labels[i].set_alpha(0.7)

    def func(label):
        index = labels.index(label)
        lines[index].set_visible(not lines[index].get_visible())
        plt.draw()    

    check.on_clicked(func)
    plt.show()

def GerarFatoresCompressibilidade(tfn_fc, lbrR_fc, tfT_fc, tfa_fc, tfb_fc, uV):
    frToolBarBottom.grid_forget()
    Rv = lbrR_fc.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])

    unidadeVolume = uV.get()
    n = float(tfn_fc.get())
    T = float(tfT_fc.get())
    a = float(tfa_fc.get())
    b = float(tfb_fc.get())

    V = np.arange(0.1, 5, 0.01)
    s0 = ((V)/(V-n*b))-((n*a)/(V*R*T))
    s1 = ((V)/(V-n*b))-((n*a)/(V*R*T**2))
    s2 = ((V)/(V-n*b))-((n*a)/((V+n*b)*R*T**(3/2)))
    s3 = ((V)/(V-n*b))*(np.exp((-n*a)/(R*T*V)))
    s4 = V/V

    fig, ax = plt.subplots()
    l0, = ax.plot(V, s0, lw=2, color='k', label='VdW')
    l1, = ax.plot(V, s1, visible=False, lw=2, color='r', label='Berth')
    l2, = ax.plot(V, s2, visible=False, lw=2, color='g', label='RK')
    l3, = ax.plot(V, s3, visible=False, lw=2, color='y', label='Diet')
    l4, = ax.plot(V, s4, visible=False, lw=2, color='b', label='Ideal')
    plt.subplots_adjust(left=0.3)

    ax.set_ylim(bottom=0., top = 20)


    plt.title("Fator de Compressibilidade Z(V)")
    plt.xlabel("V "+"("+unidadeVolume+")")
    plt.ylabel("Z - Fator de Compressibilidade")
    plt.grid(alpha = .4, linestyle = "--")

    lines = [l0, l1, l2, l3, l4]

    # Make checkbuttons with all plotted lines with correct visibility
    rax = plt.axes([0.05, 0.4, 0.1, 0.15])
    labels = [str(line.get_label()) for line in lines]
    visibility = [line.get_visible() for line in lines]
    check = CheckButtons(rax, labels, visibility)
    for i, c in enumerate(["k", "r", "g", "y", "b"]):
        check.labels[i].set_color(c)
        check.labels[i].set_alpha(0.7)

    def func(label):
        index = labels.index(label)
        lines[index].set_visible(not lines[index].get_visible())
        plt.draw()    

    check.on_clicked(func)
    plt.show()

def TratamentoDados(raizes):
    for i in range(len(raizes)):
        if raizes[i] > 0:
            if "j" in str(raizes[i]) and "+" in str(raizes[i]):
                partes = str(raizes[i]).split("+")
                if float(partes[1][0:len(partes[1])-2]) == 0:
                    return float(partes[0][1:len(partes[0])])
            elif "j" in str(raizes[i]) and "-" in str(raizes[i]):
                partes = str(raizes[i]).split("-")
                if float(partes[1][0:len(partes[1])-2]) == 0:
                    return float(partes[0][1:len(partes[0])])
            else:
                return raizes[i]

def Regula_Falsi(nomeGrandezas, grandezas, funcao):
    erro = 0.0001
    f = lambda x : eval(funcao)
    x0 = 0
    x1 = 0

    a = grandezas[4]
    b = grandezas[5]

    if not("V" in nomeGrandezas):
        P = grandezas[0]
        n = grandezas[1]
        R = grandezas[2]
        T = grandezas[3]

        x0 = (n*R*T)/(2*P)
        x1 = (3*n*R*T)/(2*P)
        
    elif not("n" in nomeGrandezas):
        P = grandezas[0]
        V = grandezas[1]
        R = grandezas[2]
        T = grandezas[3]

        x0 = (P*V)/(2*R*T)
        x1 = (3*P*V)/(2*R*T)

    elif not("T" in nomeGrandezas):
        P = grandezas[0]
        V = grandezas[1]
        n = grandezas[2]
        R = grandezas[3]

        x0 = (P*V)/(2*n*R)
        x1 = (3*P*V)/(2*n*R)

    if f(x0)*f(x1) < 0:
        while True:
            x2 = (x0*f(x1)-x1*f(x0))/(f(x1)-f(x0))
            if abs((x2-x0)/x2) < erro or abs((x2-x1)/x2) < erro:
                break
            else:
                if f(x2)*f(x0) < 0:
                    x0 = x0
                    x1 = x2
                elif f(x2)*f(x1) < 0:
                    x0 = x2
                    x1 = x1
        return x2

def EquacaoDeDietericiCalculo(tfP_d, tfV_d, tfn_d, lbrR_d, tfT_d, tfa_d, tfb_d, uP, uV):
    frToolBarBottom.grid_forget()
    UnidadesA = lba.cget("text")
    Rv = lbrR_d.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])
    unidadePressao = uP.get()
    unidadeVolume = uV.get()

    if len(tfP_d.get()) == 0 and len(tfV_d.get()) == 0 and len(tfn_d.get()) == 0 and len(tfT_d.get()) and len(tfa_d.get()) and len(tfb_d.get()) == 0:
        lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
        frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
        frToolBarBottom.config(bg = "red")
        lbAviso.grid(row = 0, column = 0)
        frToolBarBottom.columnconfigure(0, weight=1)
        frToolBarBottom.rowconfigure(0, weight = 1)
    elif len(tfP_d.get()) != 0 and len(tfV_d.get()) != 0 and len(tfn_d.get()) != 0 and len(tfT_d.get()) != 0 and len(tfa_d.get()) != 0 and len(tfb_d.get()) != 0:
        P = float(tfP_d.get())
        V = float(tfV_d.get())
        n = float(tfn_d.get())
        T = float(tfT_d.get())
        a = float(tfa_d.get())
        b = float(tfb_d.get())

        if P*(V-n*b) == n*R*T*exp((-n*a)/(R*T*V)):
            lbAviso = Label(frToolBarBottom, text = "Igualdade válida", bg = "green", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            lbAviso = Label(frToolBarBottom, text = "Igualdade inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfP_d.get()) == 0:
        if len(tfV_d.get()) == 0 or len(tfn_d.get()) == 0 or len(tfT_d.get()) == 0 or len(tfa_d.get()) == 0 or len(tfb_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            V = float(tfV_d.get())
            n = float(tfn_d.get())
            T = float(tfT_d.get())
            a = float(tfa_d.get())
            b = float(tfb_d.get())
            P = (n*R*T*exp((-n*a)/(R*T*V)))/(V-n*b)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbPDesejado = Label(frToolBarBottom, text = "P = "+str(P)+" "+unidadePressao, bg = "green", width = 40)
            lbPDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)
                
    elif len(tfV_d.get()) == 0:
        if len(tfP_d.get()) == 0 or len(tfn_d.get()) == 0 or len(tfT_d.get()) == 0 or len(tfa_d.get()) == 0 or len(tfb_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_d.get())
            n = float(tfn_d.get())
            T = float(tfT_d.get())
            a = float(tfa_d.get())
            b = float(tfb_d.get())

            #FUNÇÃO PARA O MÉTODO NUMÉRICO
            funcao = str(P)+"*x-"+str(n*R*T)+"*exp((-"+str(n*a)+")/("+str(R*T)+"*x))-"+str(P*n*b)
            nomeGrandezas = ["P", "n", "R", "T", "a", "b"]
            grandezas = [P, n, R, T, a, b]
            V = Regula_Falsi(nomeGrandezas, grandezas, funcao)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbVDesejado = Label(frToolBarBottom, text = "V = "+str(V)+" "+unidadeVolume, bg = "green", width = 40)
            lbVDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfn_d.get()) == 0:
        if len(tfP_d.get()) == 0 or len(tfV_d.get()) == 0 or len(tfT_d.get()) == 0 or len(tfa_d.get()) == 0 or len(tfb_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_d.get())
            V = float(tfV_d.get())
            T = float(tfT_d.get())
            a = float(tfa_d.get())
            b = float(tfb_d.get())
            
            #FUNÇÃO PARA O MÉTODO NUMÉRICO
            funcao = str(P*V)+"-x*"+str(R*T)+"*exp((-x*"+str(a)+")/("+str(R*T*V)+"))-x*"+str(P*b)
            nomeGrandezas = ["P", "V", "R", "T", "a", "b"]
            grandezas = [P, V, R, T, a, b]
            n = Regula_Falsi(nomeGrandezas, grandezas, funcao)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbnDesejado = Label(frToolBarBottom, text = "n = "+str(n)+" mol(s)", bg = "green", width = 40)
            lbnDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfT_d.get()) == 0:
        if len(tfP_d.get()) == 0 or len(tfV_d.get()) == 0 or len(tfn_d.get()) == 0 or len(tfa_d.get()) == 0 or len(tfb_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_d.get())
            V = float(tfV_d.get())
            n = float(tfn_d.get())
            a = float(tfa_d.get())
            b = float(tfb_d.get())

            #FUNÇÃO PARA O MÉTODO NUMÉRICO
            funcao = "P*V-n*R*x*exp((-n*a)/(R*x*V))-P*n*b"
            funcao = str(P*V)+"-x*"+str(n*R)+"*exp((-"+str(n*a)+")/("+str(R*V)+"*x))-"+str(P*n*b)
            nomeGrandezas = ["P", "V", "n", "R", "a", "b"]
            grandezas = [P, V, n, R, a, b]
            T = Regula_Falsi(nomeGrandezas, grandezas, funcao)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "T = "+str(T)+" K", bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfa_d.get()) == 0:
        if len(tfP_d.get()) == 0 or len(tfV_d.get()) == 0 or len(tfn_d.get()) == 0 or len(tfT_d.get()) == 0 or len(tfb_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_d.get())
            V = float(tfV_d.get())
            n = float(tfn_d.get())
            b = float(tfb_d.get())
            T = float(tfT_d.get())

            a = (-(R*T*V)/n)*log((P*(V-n*b))/(n*R*T))

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "a = "+str(a)+" "+UnidadesA, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfb_d.get()) == 0:
        if len(tfP_d.get()) == 0 or len(tfV_d.get()) == 0 or len(tfn_d.get()) == 0 or len(tfT_d.get()) == 0 or len(tfa_d.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_d.get())
            V = float(tfV_d.get())
            n = float(tfn_d.get())
            a = float(tfa_d.get())
            T = float(tfT_d.get())
            b = (1/n)*(V-(n*R*T*exp((-n*a)/(R*T*V)))/P)
            
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "b = "+str(b)+" "+unidadeVolume, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

def EquacaoDeRedlichKwongCalculo(tfP_rk, tfV_rk, tfn_rk, lbrR_rk, tfT_rk, tfa_rk, tfb_rk, uP, uV):
    frToolBarBottom.grid_forget()
    UnidadesA = lba.cget("text")
    Rv = lbrR_rk.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])
    unidadePressao = uP.get()
    unidadeVolume = uV.get()

    if len(tfP_rk.get()) == 0 and len(tfV_rk.get()) == 0 and len(tfn_rk.get()) == 0 and len(tfT_rk.get()) and len(tfa_rk.get()) and len(tfb_rk.get()) == 0:
        lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
        frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
        frToolBarBottom.config(bg = "red")
        lbAviso.grid(row = 0, column = 0)
        frToolBarBottom.columnconfigure(0, weight=1)
        frToolBarBottom.rowconfigure(0, weight = 1)
    elif len(tfP_rk.get()) != 0 and len(tfV_rk.get()) != 0 and len(tfn_rk.get()) != 0 and len(tfT_rk.get()) != 0 and len(tfa_rk.get()) != 0 and len(tfb_rk.get()) != 0:
        P = float(tfP_rk.get())
        V = float(tfV_rk.get())
        n = float(tfn_rk.get())
        T = float(tfT_rk.get())
        a = float(tfa_rk.get())
        b = float(tfb_rk.get())

        if (P+(a*n**2)/((T**(1/2))*V*(V+n*b)))*(V-n*b) == n*R*T:
            lbAviso = Label(frToolBarBottom, text = "Igualdade válida", bg = "green", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            lbAviso = Label(frToolBarBottom, text = "Igualdade inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfP_rk.get()) == 0:
            if len(tfV_rk.get()) == 0 or len(tfn_rk.get()) == 0 or len(tfT_rk.get()) == 0 or len(tfa_rk.get()) == 0 or len(tfb_rk.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                V = float(tfV_rk.get())
                n = float(tfn_rk.get())
                T = float(tfT_rk.get())
                a = float(tfa_rk.get())
                b = float(tfb_rk.get())
                P = (n*R*T)/(V-n*b) - (a*n**2)/((T**(1/2))*V*(V+n*b))

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbPDesejado = Label(frToolBarBottom, text = "P = "+str(P)+" "+unidadePressao, bg = "green", width = 40)
                lbPDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)
                
    elif len(tfV_rk.get()) == 0:
            if len(tfP_rk.get()) == 0 or len(tfn_rk.get()) == 0 or len(tfT_rk.get()) == 0 or len(tfa_rk.get()) == 0 or len(tfb_rk.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_rk.get())
                n = float(tfn_rk.get())
                T = float(tfT_rk.get())
                a = float(tfa_rk.get())
                b = float(tfb_rk.get())

                #Coeficientes do polinômio
                A = P*T**(1/2)
                B = -(n*R*T**(3/2))
                C = a*n**2 - P*n**2*b**2*T**(1/2)-n**2*b*R*T**(3/2)
                D = -a*b*n**3
                r = np.roots([A, B, C, D])
                V = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbVDesejado = Label(frToolBarBottom, text = "V = "+str(V)+" "+unidadeVolume, bg = "green", width = 40)
                lbVDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfn_rk.get()) == 0:
            if len(tfP_rk.get()) == 0 or len(tfV_rk.get()) == 0 or len(tfT_rk.get()) == 0 or len(tfa_rk.get()) == 0 or len(tfb_rk.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_rk.get())
                V = float(tfV_rk.get())
                T = float(tfT_rk.get())
                a = float(tfa_rk.get())
                b = float(tfb_rk.get())
                
                #Coeficientes do polinômio
                A = a*b
                B = R*b*V*T**(3/2) + P*V*b**2*T**(1/2) - a*V
                C = R*V**2*T**(3/2)
                D = -P*V**3*T**(1/2)
                r = np.roots([A, B, C, D])
                n = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbnDesejado = Label(frToolBarBottom, text = "n = "+str(n)+" mol(s)", bg = "green", width = 40)
                lbnDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfT_rk.get()) == 0:
            if len(tfP_rk.get()) == 0 or len(tfV_rk.get()) == 0 or len(tfn_rk.get()) == 0 or len(tfa_rk.get()) == 0 or len(tfb_rk.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_rk.get())
                V = float(tfV_rk.get())
                n = float(tfn_rk.get())
                a = float(tfa_rk.get())
                b = float(tfb_rk.get())

                #COEFICIENTES DO POLINOMIO
                A = (n*R*V**2+n**2*R*b*V)**2
                B = 2*(n*R*V**2+n**2*R*b*V)*(P*V*n**2*b**2-P*V**3)
                C = (P*V*n**2*b**2-P*V**3)**2
                D = -(a*V*n**2-a*b*n**3)**2
                r = np.roots([A, B, C, D])
                print(r)
                T = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbTDesejado = Label(frToolBarBottom, text = "T = "+str(T)+" K", bg = "green", width = 40)
                lbTDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfa_rk.get()) == 0:
        if len(tfP_rk.get()) == 0 or len(tfV_rk.get()) == 0 or len(tfn_rk.get()) == 0 or len(tfT_rk.get()) == 0 or len(tfb_rk.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_rk.get())
            V = float(tfV_rk.get())
            n = float(tfn_rk.get())
            b = float(tfb_rk.get())
            T = float(tfT_rk.get())

            a = (((T**(1/2))*V*(V+n*b))/(n**2))*(((n*R*T)/(V-n*b))-P)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "a = "+str(a)+" "+UnidadesA, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfb_rk.get()) == 0:
        if len(tfP_rk.get()) == 0 or len(tfV_rk.get()) == 0 or len(tfn_rk.get()) == 0 or len(tfT_rk.get()) == 0 or len(tfa_rk.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_rk.get())
            V = float(tfV_rk.get())
            n = float(tfn_rk.get())
            a = float(tfa_rk.get())
            T = float(tfT_rk.get())

            #COEFICIENTES DO POLINOMIO
            A = P*V*n**2*T**(1/2)
            B = R*V*n**2*T**(3/2) + a*n**3
            C = n*R*V**2*T**(3/2)-P*V**3*T**(1/2)-a*V*n**2
            r = np.roots([A, B, C])
            b = TratamentoDados(r)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "b = "+str(b)+" "+unidadeVolume, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

def EquacaoDeBerthelotCalculo(tfP_b, tfV_b, tfn_b, lbrR_b, tfT_b, tfa_b, tfb_b, uP, uV):
    frToolBarBottom.grid_forget()
    UnidadesA = lba.cget("text")
    Rv = lbrR_b.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])
    unidadePressao = uP.get()
    unidadeVolume = uV.get()

    if len(tfP_b.get()) == 0 and len(tfV_b.get()) == 0 and len(tfn_b.get()) == 0 and len(tfT_b.get()) and len(tfa_b.get()) and len(tfb_b.get()) == 0:
        lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
        frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
        frToolBarBottom.config(bg = "red")
        lbAviso.grid(row = 0, column = 0)
        frToolBarBottom.columnconfigure(0, weight=1)
        frToolBarBottom.rowconfigure(0, weight = 1)
    elif len(tfP_b.get()) != 0 and len(tfV_b.get()) != 0 and len(tfn_b.get()) != 0 and len(tfT_b.get()) != 0 and len(tfa_b.get()) != 0 and len(tfb_b.get()) != 0:
        P = float(tfP_b.get())
        V = float(tfV_b.get())
        n = float(tfn_b.get())
        T = float(tfT_b.get())
        a = float(tfa_b.get())
        b = float(tfb_b.get())

        if (P+(a*n**2)/(T*V**2))*(V-n*b) == n*R*T:
            lbAviso = Label(frToolBarBottom, text = "Igualdade válida", bg = "green", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            lbAviso = Label(frToolBarBottom, text = "Igualdade inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfP_b.get()) == 0:
            if len(tfV_b.get()) == 0 or len(tfn_b.get()) == 0 or len(tfT_b.get()) == 0 or len(tfa_b.get()) == 0 or len(tfb_b.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                V = float(tfV_b.get())
                n = float(tfn_b.get())
                T = float(tfT_b.get())
                a = float(tfa_b.get())
                b = float(tfb_b.get())
                P = (n*R*T)/(V-n*b) - (a*n**2)/(T*V**2)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbPDesejado = Label(frToolBarBottom, text = "P = "+str(P)+" "+unidadePressao, bg = "green", width = 40)
                lbPDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)
                
    elif len(tfV_b.get()) == 0:
            if len(tfP_b.get()) == 0 or len(tfn_b.get()) == 0 or len(tfT_b.get()) == 0 or len(tfa_b.get()) == 0 or len(tfb_b.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_b.get())
                n = float(tfn_b.get())
                T = float(tfT_b.get())
                a = float(tfa_b.get())
                b = float(tfb_b.get())

                #Coeficientes do polinômio
                A = P*T
                B = -(P*T*n*b+n*R*T**2)
                C = a*n**2
                D = -a*b*n**3
                r = np.roots([A, B, C, D])
                V = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbVDesejado = Label(frToolBarBottom, text = "V = "+str(V)+" "+unidadeVolume, bg = "green", width = 40)
                lbVDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfn_b.get()) == 0:
            if len(tfP_b.get()) == 0 or len(tfV_b.get()) == 0 or len(tfT_b.get()) == 0 or len(tfa_b.get()) == 0 or len(tfb_b.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_b.get())
                V = float(tfV_b.get())
                T = float(tfT_b.get())
                a = float(tfa_b.get())
                b = float(tfb_b.get())
                
                #Coeficientes do polinômio
                A = a*b
                B = -a*V
                C = P*T*b*V**2+R*T**2*V**2
                D = -T*P*V**3
                r = np.roots([A, B, C, D])
                n = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbnDesejado = Label(frToolBarBottom, text = "n = "+str(n)+" mol(s)", bg = "green", width = 40)
                lbnDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfT_b.get()) == 0:
            if len(tfP_b.get()) == 0 or len(tfV_b.get()) == 0 or len(tfn_b.get()) == 0 or len(tfa_b.get()) == 0 or len(tfb_b.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_b.get())
                V = float(tfV_b.get())
                n = float(tfn_b.get())
                a = float(tfa_b.get())
                b = float(tfb_b.get())

                #COEFICIENTES DO POLINOMIO
                A = n*R*V**2
                B = n*b*P*V**2 - P*V**3
                C = a*b*n**3 - a*V*n**2
                r = np.roots([A, B, C])
                T = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbTDesejado = Label(frToolBarBottom, text = "T = "+str(T)+" K", bg = "green", width = 40)
                lbTDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfa_b.get()) == 0:
        if len(tfP_b.get()) == 0 or len(tfV_b.get()) == 0 or len(tfn_b.get()) == 0 or len(tfT_b.get()) == 0 or len(tfb_b.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_b.get())
            V = float(tfV_b.get())
            n = float(tfn_b.get())
            b = float(tfb_b.get())
            T = float(tfT_b.get())

            a = ((T*V**2)/n**2)*((n*R*T)/(V-n*b)-P)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "a = "+str(a)+" "+UnidadesA, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfb_b.get()) == 0:
        if len(tfP_b.get()) == 0 or len(tfV_b.get()) == 0 or len(tfn_b.get()) == 0 or len(tfT_b.get()) == 0 or len(tfa_b.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_b.get())
            V = float(tfV_b.get())
            n = float(tfn_b.get())
            a = float(tfa_b.get())
            T = float(tfT_b.get())
            b = (1/n)*(V-(n*R*T)/(P+(a*n**2)/(T*V**2)))
            
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "b = "+str(b)+" "+unidadeVolume, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)
            
def EquacaoDeVanDerWaalsCalculo(tfP_vdw, tfV_vdw, tfn_vdw, lbrR, tfT_vdw, tfa_vdw, tfb_vdw, uP, uV):
    frToolBarBottom.grid_forget()
    UnidadesA = lba.cget("text")
    Rv = lbrR.cget("text").split(" ")
    if "x10^" in Rv[0]:
        valores = Rv[0].split("x10^")
        final = valores[0]+"e"+valores[1]
        R = float(final)
    else:
        R = float(Rv[0])
    unidadePressao = uP.get()
    unidadeVolume = uV.get()

    if len(tfP_vdw.get()) == 0 and len(tfV_vdw.get()) == 0 and len(tfn_vdw.get()) == 0 and len(tfT_vdw.get()) and len(tfa_vdw.get()) and len(tfb_vdw.get()) == 0:
        lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
        frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
        frToolBarBottom.config(bg = "red")
        lbAviso.grid(row = 0, column = 0)
        frToolBarBottom.columnconfigure(0, weight=1)
        frToolBarBottom.rowconfigure(0, weight = 1)
    elif len(tfP_vdw.get()) != 0 and len(tfV_vdw.get()) != 0 and len(tfn_vdw.get()) != 0 and len(tfT_vdw.get()) != 0 and len(tfa_vdw.get()) != 0 and len(tfb_vdw.get()) != 0:
        P = float(tfP_vdw.get())
        V = float(tfV_vdw.get())
        n = float(tfn_vdw.get())
        T = float(tfT_vdw.get())
        a = float(tfa_vdw.get())
        b = float(tfb_vdw.get())

        if (P+(a*n**2)/V**2)*(V-n*b) == n*R*T:
            lbAviso = Label(frToolBarBottom, text = "Igualdade válida", bg = "green", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            lbAviso = Label(frToolBarBottom, text = "Igualdade inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfP_vdw.get()) == 0:
            if len(tfV_vdw.get()) == 0 or len(tfn_vdw.get()) == 0 or len(tfT_vdw.get()) == 0 or len(tfa_vdw.get()) == 0 or len(tfb_vdw.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                V = float(tfV_vdw.get())
                n = float(tfn_vdw.get())
                T = float(tfT_vdw.get())
                a = float(tfa_vdw.get())
                b = float(tfb_vdw.get())
                P = (n*R*T)/(V-n*b) - (a*n**2)/V**2

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbPDesejado = Label(frToolBarBottom, text = "P = "+str(P)+" "+unidadePressao, bg = "green", width = 40)
                lbPDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)
                
    elif len(tfV_vdw.get()) == 0:
            if len(tfP_vdw.get()) == 0 or len(tfn_vdw.get()) == 0 or len(tfT_vdw.get()) == 0 or len(tfa_vdw.get()) == 0 or len(tfb_vdw.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_vdw.get())
                n = float(tfn_vdw.get())
                T = float(tfT_vdw.get())
                a = float(tfa_vdw.get())
                b = float(tfb_vdw.get())

                #Coeficientes do polinômio
                A = P
                B = -(n*b*P+n*R*T)
                C = a*n**2
                D = -a*b*n**3
                r = np.roots([A, B, C, D])
                V = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbVDesejado = Label(frToolBarBottom, text = "V = "+str(V)+" "+unidadeVolume, bg = "green", width = 40)
                lbVDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfn_vdw.get()) == 0:
            if len(tfP_vdw.get()) == 0 or len(tfV_vdw.get()) == 0 or len(tfT_vdw.get()) == 0 or len(tfa_vdw.get()) == 0 or len(tfb_vdw.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_vdw.get())
                V = float(tfV_vdw.get())
                T = float(tfT_vdw.get())
                a = float(tfa_vdw.get())
                b = float(tfb_vdw.get())
                
                #Coeficientes do polinômio
                A = a*b
                B = -a*V
                C = P*b*V**2+R*T*V**2
                D = -P*V**3
                r = np.roots([A, B, C, D])
                n = TratamentoDados(r)

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbnDesejado = Label(frToolBarBottom, text = "n = "+str(n)+" mol(s)", bg = "green", width = 40)
                lbnDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfT_vdw.get()) == 0:
            if len(tfP_vdw.get()) == 0 or len(tfV_vdw.get()) == 0 or len(tfn_vdw.get()) == 0 or len(tfa_vdw.get()) == 0 or len(tfb_vdw.get()) == 0:
                lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "red")
                lbAviso.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight=1)
                frToolBarBottom.rowconfigure(0, weight = 1)
            else:
                P = float(tfP_vdw.get())
                V = float(tfV_vdw.get())
                n = float(tfn_vdw.get())
                a = float(tfa_vdw.get())
                b = float(tfb_vdw.get())
                T = (P+(a*n**2)/V**2)*((V-n*b)/(n*R))

                frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
                frToolBarBottom.config(bg = "green")
                lbTDesejado = Label(frToolBarBottom, text = "T = "+str(T)+" K", bg = "green", width = 40)
                lbTDesejado.grid(row = 0, column = 0)
                frToolBarBottom.columnconfigure(0, weight = 1)
                frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfa_vdw.get()) == 0:
        if len(tfP_vdw.get()) == 0 or len(tfV_vdw.get()) == 0 or len(tfn_vdw.get()) == 0 or len(tfT_vdw.get()) == 0 or len(tfb_vdw.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_vdw.get())
            V = float(tfV_vdw.get())
            n = float(tfn_vdw.get())
            b = float(tfb_vdw.get())
            T = float(tfT_vdw.get())

            a = (V**2/n**2)*(((n*R*T)/(V-n*b))-P)

            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "a = "+str(a)+" "+UnidadesA, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

    elif len(tfb_vdw.get()) == 0:
        if len(tfP_vdw.get()) == 0 or len(tfV_vdw.get()) == 0 or len(tfn_vdw.get()) == 0 or len(tfT_vdw.get()) == 0 or len(tfa_vdw.get()) == 0:
            lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red", width = 40)
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "red")
            lbAviso.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight=1)
            frToolBarBottom.rowconfigure(0, weight = 1)
        else:
            P = float(tfP_vdw.get())
            V = float(tfV_vdw.get())
            n = float(tfn_vdw.get())
            a = float(tfa_vdw.get())
            T = float(tfT_vdw.get())
            b = (1/n)*(V-(n*R*T)/(P+(a*n**2)/V**2))

            
            frToolBarBottom.grid(row = 2, column = 0, sticky = W+E)
            frToolBarBottom.config(bg = "green")
            lbTDesejado = Label(frToolBarBottom, text = "b = "+str(b)+" "+unidadeVolume, bg = "green", width = 40)
            lbTDesejado.grid(row = 0, column = 0)
            frToolBarBottom.columnconfigure(0, weight = 1)
            frToolBarBottom.rowconfigure(0, weight = 1)

def escolhaTela(tipo):
    lbAviso.grid_forget()
    frToolBarBottom.grid_forget()
    if vetor[0] == "VanDerWaals":
        frVanDerWaals.grid_forget()
    elif vetor[0] == "Berthelot":
        frBerthelot.grid_forget()
    elif vetor[0] == "Dieterici":
        frDieterici.grid_forget()
    elif vetor[0] == "RedlichKwong":
        frRedlichKwong.grid_forget()
    elif vetor[0] == "FatoresCompressibilidade":
        frFatorCompressibilidade.grid_forget()
    elif vetor[0] == "Isotermas":
        frIsotermas.grid_forget()

    if tipo == "van der Waals":
        vetor[0] = "VanDerWaals"
        frVanDerWaals.grid(row = 1, column = 0)
    elif tipo == "Berthelot":
        vetor[0] = "Berthelot"
        frBerthelot.grid(row = 1, column = 0)
    elif tipo == "Dieterici":
        vetor[0] = "Dieterici"
        frDieterici.grid(row = 1, column = 0)
    elif tipo == "Redlich-Kwong":
        vetor[0] = "RedlichKwong"
        frRedlichKwong.grid(row = 1, column = 0)
    elif tipo == "FatoresCompressibilidade":
        vetor[0] = "FatoresCompressibilidade"
        frFatorCompressibilidade.grid(row = 1, column = 0)
    elif tipo == "Isotermas":
        vetor[0] = "Isotermas"
        frIsotermas.grid(row = 1, column = 0)

def mudarConstante(parametro):
    if parametro == "Pa" or parametro == "atm" or parametro == "mmHg":
        if parametro == "Pa":
            aux_mudar_constante[0] = "Pa" 
            if aux_mudar_constante[1] == "m^3":
                lbrR["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31 Pa.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.Pa"
                lba_b["text"] = "m^6.Pa"
                lba_d["text"] = "m^6.Pa"
                lba_rk["text"] = "m^6.Pa"
                lba_fc["text"] = "m^6.Pa"
                lba_isot["text"] = "m^6.Pa"

            elif aux_mudar_constante[1] == "L":
                lbrR["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_b["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_d["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.Pa"
                lba_b["text"] = "L^2.Pa"
                lba_d["text"] = "L^2.Pa"
                lba_rk["text"] = "L^2.Pa"
                lba_fc["text"] = "L^2.Pa"
                lba_isot["text"] = "L^2.Pa"

            elif aux_mudar_constante[1] == "mL":
                lbrR["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_b["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_d["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.Pa"
                lba_b["text"] = "mL^2.Pa"
                lba_d["text"] = "mL^2.Pa"
                lba_rk["text"] = "mL^2.Pa"
                lba_fc["text"] = "mL^2.Pa"
                lba_isot["text"] = "mL^2.Pa"

        elif parametro == "atm":
            aux_mudar_constante[0] = "atm" 
            if aux_mudar_constante[1] == "m^3":
                lbrR["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.atm"
                lba_b["text"] = "m^6.atm"
                lba_d["text"] = "m^6.atm"
                lba_rk["text"] = "m^6.atm"
                lba_fc["text"] = "m^6.atm"
                lba_isot["text"] = "m^6.atm"

            elif aux_mudar_constante[1] == "L":
                lbrR["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_b["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_d["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_rk["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_fc["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_isot["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.atm"
                lba_b["text"] = "L^2.atm"
                lba_d["text"] = "L^2.atm"
                lba_rk["text"] = "L^2.atm"
                lba_fc["text"] = "L^2.atm"
                lba_isot["text"] = "L^2.atm"

            elif aux_mudar_constante[1] == "mL":
                lbrR["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_b["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_d["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "82 atm.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.atm"
                lba_b["text"] = "mL^2.atm"
                lba_d["text"] = "mL^2.atm"
                lba_rk["text"] = "mL^2.atm"
                lba_fc["text"] = "mL^2.atm"
                lba_isot["text"] = "mL^2.atm"

        elif parametro == "mmHg":
            aux_mudar_constante[0] = "mmHg" 
            if aux_mudar_constante[1] == "m^3":
                lbrR["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.mmHg"
                lba_b["text"] = "m^6.mmHg"
                lba_d["text"] = "m^6.mmHg"
                lba_rk["text"] = "m^6.mmHg"
                lba_fc["text"] = "m^6.mmHg"
                lba_isot["text"] = "m^6.mmHg"

            elif aux_mudar_constante[1] == "L":
                lbrR["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_b["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_d["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_rk["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_fc["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_isot["text"] = "62.32 mmHg.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.mmHg"
                lba_b["text"] = "L^2.mmHg"
                lba_d["text"] = "L^2.mmHg"
                lba_rk["text"] = "L^2.mmHg"
                lba_fc["text"] = "L^2.mmHg"
                lba_isot["text"] = "L^2.mmHg"

            elif aux_mudar_constante[1] == "mL":
                lbrR["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_b["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_d["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "62320 mmHg.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.mmHg"
                lba_b["text"] = "mL^2.mmHg"
                lba_d["text"] = "mL^2.mmHg"
                lba_rk["text"] = "mL^2.mmHg"
                lba_fc["text"] = "mL^2.mmHg"
                lba_isot["text"] = "mL^2.mmHg"

    elif parametro == "m^3" or parametro == "L" or parametro == "mL":
        if parametro == "m^3":
            aux_mudar_constante[1] = "m^3" 
            if aux_mudar_constante[0] == "Pa":
                lbrR["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31 Pa.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31 Pa.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.Pa"
                lba_b["text"] = "m^6.Pa"
                lba_d["text"] = "m^6.Pa"
                lba_rk["text"] = "m^6.Pa"
                lba_fc["text"] = "m^6.Pa"
                lba_isot["text"] = "m^6.Pa"

            elif aux_mudar_constante[0] == "atm":
                lbrR["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "8.2x10^-5 atm.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.atm"
                lba_b["text"] = "m^6.atm"
                lba_d["text"] = "m^6.atm"
                lba_rk["text"] = "m^6.atm"
                lba_fc["text"] = "m^6.atm"
                lba_isot["text"] = "m^6.atm"

            elif aux_mudar_constante[0] == "mmHg":
                lbrR["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_b["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_d["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_rk["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_fc["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"
                lbrR_isot["text"] = "6.2x10^-2 mmHg.m^3.mol^-1.K^-1"

                lbb["text"] = "m^3.mol^-1"
                lbb_b["text"] = "m^3.mol^-1"
                lbb_d["text"] = "m^3.mol^-1"
                lbb_rk["text"] = "m^3.mol^-1"
                lbb_fc["text"] = "m^3.mol^-1"
                lbb_isot["text"] = "m^3.mol^-1"

                lba["text"] = "m^6.mmHg"
                lba_b["text"] = "m^6.mmHg"
                lba_d["text"] = "m^6.mmHg"
                lba_rk["text"] = "m^6.mmHg"
                lba_fc["text"] = "m^6.mmHg"
                lba_isot["text"] = "m^6.mmHg"

        elif parametro == "L":
            aux_mudar_constante[1] = "L" 
            if aux_mudar_constante[0] == "Pa":
                lbrR["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_b["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_d["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31x10^3 Pa.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.Pa"
                lba_b["text"] = "L^2.Pa"
                lba_d["text"] = "L^2.Pa"
                lba_rk["text"] = "L^2.Pa"
                lba_fc["text"] = "L^2.Pa"
                lba_isot["text"] = "L^2.Pa"

            if aux_mudar_constante[0] == "atm":
                lbrR["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_b["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_d["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_rk["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_fc["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"
                lbrR_isot["text"] = "8.2x10^-2 atm.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.atm"
                lba_b["text"] = "L^2.atm"
                lba_d["text"] = "L^2.atm"
                lba_rk["text"] = "L^2.atm"
                lba_fc["text"] = "L^2.atm"
                lba_isot["text"] = "L^2.atm"

            if aux_mudar_constante[0] == "mmHg":
                lbrR["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_b["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_d["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_rk["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_fc["text"] = "62.32 mmHg.L.mol^-1.K^-1"
                lbrR_isot["text"] = "62.32 mmHg.L.mol^-1.K^-1"

                lbb["text"] = "L.mol^-1"
                lbb_b["text"] = "L.mol^-1"
                lbb_d["text"] = "L.mol^-1"
                lbb_rk["text"] = "L.mol^-1"
                lbb_fc["text"] = "L.mol^-1"
                lbb_isot["text"] = "L.mol^-1"

                lba["text"] = "L^2.mmHg"
                lba_b["text"] = "L^2.mmHg"
                lba_d["text"] = "L^2.mmHg"
                lba_rk["text"] = "L^2.mmHg"
                lba_fc["text"] = "L^2.mmHg"
                lba_isot["text"] = "L^2.mmHg"

        elif parametro == "mL":
            aux_mudar_constante[1] = "mL" 
            if aux_mudar_constante[0] == "Pa":
                lbrR["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_b["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_d["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "8.31x10^6 Pa.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.Pa"
                lba_b["text"] = "mL^2.Pa"
                lba_d["text"] = "mL^2.Pa"
                lba_rk["text"] = "mL^2.Pa"
                lba_fc["text"] = "mL^2.Pa"
                lba_isot["text"] = "mL^2.Pa"

            if aux_mudar_constante[0] == "atm":
                lbrR["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_b["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_d["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "82 atm.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "82 atm.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.atm"
                lba_b["text"] = "mL^2.atm"
                lba_d["text"] = "mL^2.atm"
                lba_rk["text"] = "mL^2.atm"
                lba_fc["text"] = "mL^2.atm"
                lba_isot["text"] = "mL^2.atm"

            if aux_mudar_constante[0] == "mmHg":
                lbrR["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_b["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_d["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_rk["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_fc["text"] = "62320 mmHg.mL.mol^-1.K^-1"
                lbrR_isot["text"] = "62320 mmHg.mL.mol^-1.K^-1"

                lbb["text"] = "mL.mol^-1"
                lbb_b["text"] = "mL.mol^-1"
                lbb_d["text"] = "mL.mol^-1"
                lbb_rk["text"] = "mL.mol^-1"
                lbb_fc["text"] = "mL.mol^-1"
                lbb_isot["text"] = "mL.mol^-1"

                lba["text"] = "mL^2.mmHg"
                lba_b["text"] = "mL^2.mmHg"
                lba_d["text"] = "mL^2.mmHg"
                lba_rk["text"] = "mL^2.mmHg"
                lba_fc["text"] = "mL^2.mmHg"
                lba_isot["text"] = "mL^2.mmHg"

janela = Tk()

#define the objects of REDLICH-KWONG equation
frRedlichKwong = Frame(janela)

lbVazio = Label(frRedlichKwong)
lbVazio1 = Label(frRedlichKwong)
lbVazio2 = Label(frRedlichKwong)

lbP_rk = Label(frRedlichKwong, text = "P")
tfP_rk = Entry(frRedlichKwong)

lbV_rk = Label(frRedlichKwong, text = "V")
tfV_rk = Entry(frRedlichKwong)

lbn_rk = Label(frRedlichKwong, text = "n")
tfn_rk = Entry(frRedlichKwong)

lbR_rk = Label(frRedlichKwong, text = "R")

lbT_rk = Label(frRedlichKwong, text = "T")
tfT_rk = Entry(frRedlichKwong)

lbA_rk = Label(frRedlichKwong, text = "a")
tfA_rk = Entry(frRedlichKwong)

lbB_rk = Label(frRedlichKwong, text = "b")
tfB_rk = Entry(frRedlichKwong)

tkP_rk = StringVar(frRedlichKwong)
escolhasP_rk = {"Pa", "atm", "mmHg"}
tkP_rk.set("Pa")
    
tkV_rk = StringVar(frRedlichKwong)
escolhasV_rk = {"m^3", "L", "mL"}
tkV_rk.set("m^3")

tkn_rk = StringVar(frRedlichKwong)
escolhasrk = {"mol"}
tkn_rk.set("mol")
    
tkT_rk = StringVar(frRedlichKwong)
escolhasT_rk = {"K"}
tkT_rk.set("K")

lbUnidades_rk = Label(frRedlichKwong, text = "Selecione as unidades da entrada")

lbPU_rk = Label(frRedlichKwong, text = "Pressão")
mnP_rk = OptionMenu(frRedlichKwong, tkP_rk, *escolhasP_rk, command = mudarConstante)
mnP_rk.config(width = 18)

lbVU_rk = Label(frRedlichKwong, text = "Volume")
mnV_rk = OptionMenu(frRedlichKwong, tkV_rk, *escolhasV_rk, command = mudarConstante)
mnV_rk.config(width = 18)

lbnU_rk = Label(frRedlichKwong, text = "Mols")
mnn_rk = OptionMenu(frRedlichKwong, tkn_rk, *escolhasrk)
mnn_rk.config(width = 18)

lbTU_rk = Label(frRedlichKwong, text = "Temperatura")
mnT_rk = OptionMenu(frRedlichKwong, tkT_rk, *escolhasT_rk)
mnT_rk.config(width = 18)

lbrR_rk = Label(frRedlichKwong, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frRedlichKwong, text = "a")
lba_rk = Label(frRedlichKwong, text = "m^6.Pa")
lbbb = Label(frRedlichKwong, text = "b")
lbb_rk = Label(frRedlichKwong, text = "m^3.mol^-1")

btCalcularRedlichKwong = Button(frRedlichKwong, text = "Calcular", width = 18)
btCalcularRedlichKwong["command"] = partial(EquacaoDeRedlichKwongCalculo,tfP_rk, tfV_rk, tfn_rk, lbrR_rk, tfT_rk, tfA_rk, tfB_rk, tkP_rk, tkV_rk)

#define the positions of the objects in the REDLICH-KWONG
lbVazio2.grid(row = 0, column = 0)
    
lbP_rk.grid(row = 1, column = 0)
tfP_rk.grid(row = 1, column = 1)
lbV_rk.grid(row = 2, column = 0)
tfV_rk.grid(row = 2, column = 1)
lbn_rk.grid(row = 3, column = 0)
tfn_rk.grid(row = 3, column = 1)
lbR_rk.grid(row = 4, column = 0)
lbrR_rk.grid(row = 4, column = 1)
lbT_rk.grid(row = 5, column = 0)
tfT_rk.grid(row = 5, column = 1)
lbA_rk.grid(row = 6, column = 0)
tfA_rk.grid(row = 6, column = 1)
lbB_rk.grid(row = 7, column = 0)
tfB_rk.grid(row = 7, column = 1)

lbVazio1.grid(row = 8, column = 0)

btCalcularRedlichKwong.grid(row = 9, column = 1)

lbUnidades_rk.grid(row = 10, column = 0)
lbPU_rk.grid(row = 11, column = 0)
mnP_rk.grid(row = 11, column = 1)
lbVU_rk.grid(row = 12, column = 0)
mnV_rk.grid(row = 12, column = 1)
lbnU_rk.grid(row = 13, column = 0)
mnn_rk.grid(row = 13, column = 1)
lbR_rk.grid(row = 14, column = 0)
lbrR_rk.grid(row = 14, column = 1)
lbba.grid(row = 15, column = 0)
lba_rk.grid(row = 15, column = 1)
lbbb.grid(row = 16, column = 0)
lbb_rk.grid(row = 16, column = 1)
lbTU_rk.grid(row = 17, column = 0)
mnT_rk.grid(row = 17, column = 1)


#define the objects of DIETERICI equation
frDieterici = Frame(janela)

lbVazio = Label(frDieterici)
lbVazio1 = Label(frDieterici)
lbVazio2 = Label(frDieterici)

lbP_d = Label(frDieterici, text = "P")
tfP_d = Entry(frDieterici)

lbV_d = Label(frDieterici, text = "V")
tfV_d = Entry(frDieterici)

lbn_d = Label(frDieterici, text = "n")
tfn_d = Entry(frDieterici)

lbR_d = Label(frDieterici, text = "R")

lbT_d = Label(frDieterici, text = "T")
tfT_d = Entry(frDieterici)

lbA_d = Label(frDieterici, text = "a")
tfA_d = Entry(frDieterici)

lbB_d = Label(frDieterici, text = "b")
tfB_d = Entry(frDieterici)

tkP_d = StringVar(frDieterici)
escolhasP_d = {"Pa", "atm", "mmHg"}
tkP_d.set("Pa")
    
tkV_d = StringVar(frDieterici)
escolhasV_d = {"m^3", "L", "mL"}
tkV_d.set("m^3")

tkn_d = StringVar(frDieterici)
escolhasd = {"mol"}
tkn_d.set("mol")
    
tkT_d = StringVar(frDieterici)
escolhasT_d = {"K"}
tkT_d.set("K")

lbUnidades_d = Label(frDieterici, text = "Selecione as unidades da entrada")

lbPU_d = Label(frDieterici, text = "Pressão")
mnP_d = OptionMenu(frDieterici, tkP_d, *escolhasP_d, command = mudarConstante)
mnP_d.config(width = 18)

lbVU_d = Label(frDieterici, text = "Volume")
mnV_d = OptionMenu(frDieterici, tkV_d, *escolhasV_d, command = mudarConstante)
mnV_d.config(width = 18)

lbnU_d = Label(frDieterici, text = "Mols")
mnn_d = OptionMenu(frDieterici, tkn_d, *escolhasd)
mnn_d.config(width = 18)

lbTU_d = Label(frDieterici, text = "Temperatura")
mnT_d = OptionMenu(frDieterici, tkT_d, *escolhasT_d)
mnT_d.config(width = 18)

lbrR_d = Label(frDieterici, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frDieterici, text = "a")
lba_d = Label(frDieterici, text = "m^6.Pa")
lbbb = Label(frDieterici, text = "b")
lbb_d = Label(frDieterici, text = "m^3.mol^-1")

btCalcularDieterici = Button(frDieterici, text = "Calcular", width = 18)
btCalcularDieterici["command"] = partial(EquacaoDeDietericiCalculo,tfP_d, tfV_d, tfn_d, lbrR_d, tfT_d, tfA_d, tfB_d, tkP_d, tkV_d)

#define the positions of the objects in the DIETERICI
lbVazio2.grid(row = 0, column = 0)
    
lbP_d.grid(row = 1, column = 0)
tfP_d.grid(row = 1, column = 1)
lbV_d.grid(row = 2, column = 0)
tfV_d.grid(row = 2, column = 1)
lbn_d.grid(row = 3, column = 0)
tfn_d.grid(row = 3, column = 1)
lbR_d.grid(row = 4, column = 0)
lbrR_d.grid(row = 4, column = 1)
lbT_d.grid(row = 5, column = 0)
tfT_d.grid(row = 5, column = 1)
lbA_d.grid(row = 6, column = 0)
tfA_d.grid(row = 6, column = 1)
lbB_d.grid(row = 7, column = 0)
tfB_d.grid(row = 7, column = 1)

lbVazio1.grid(row = 8, column = 0)

btCalcularDieterici.grid(row = 9, column = 1)

lbUnidades_d.grid(row = 10, column = 0)
lbPU_d.grid(row = 11, column = 0)
mnP_d.grid(row = 11, column = 1)
lbVU_d.grid(row = 12, column = 0)
mnV_d.grid(row = 12, column = 1)
lbnU_d.grid(row = 13, column = 0)
mnn_d.grid(row = 13, column = 1)
lbR_d.grid(row = 14, column = 0)
lbrR_d.grid(row = 14, column = 1)
lbba.grid(row = 15, column = 0)
lba_d.grid(row = 15, column = 1)
lbbb.grid(row = 16, column = 0)
lbb_d.grid(row = 16, column = 1)
lbTU_d.grid(row = 17, column = 0)
mnT_d.grid(row = 17, column = 1)

#define the objects of BERTHELOT equation
frBerthelot = Frame(janela)

lbVazio = Label(frBerthelot)
lbVazio1 = Label(frBerthelot)
lbVazio2 = Label(frBerthelot)

lbP_b = Label(frBerthelot, text = "P")
tfP_b = Entry(frBerthelot)

lbV_b = Label(frBerthelot, text = "V")
tfV_b = Entry(frBerthelot)

lbn_b = Label(frBerthelot, text = "n")
tfn_b = Entry(frBerthelot)

lbR_b = Label(frBerthelot, text = "R")

lbT_b = Label(frBerthelot, text = "T")
tfT_b = Entry(frBerthelot)

lbA_b = Label(frBerthelot, text = "a")
tfA_b = Entry(frBerthelot)

lbB_b = Label(frBerthelot, text = "b")
tfB_b = Entry(frBerthelot)

tkP_b = StringVar(frBerthelot)
escolhasP_b = {"Pa", "atm", "mmHg"}
tkP_b.set("Pa")
    
tkV_b = StringVar(frBerthelot)
escolhasV_b = {"m^3", "L", "mL"}
tkV_b.set("m^3")

tkn_b = StringVar(frBerthelot)
escolhasb = {"mol"}
tkn_b.set("mol")
    
tkT_b = StringVar(frBerthelot)
escolhasT_b = {"K"}
tkT_b.set("K")

lbUnidades_b = Label(frBerthelot, text = "Selecione as unidades da entrada")

lbPU_b = Label(frBerthelot, text = "Pressão")
mnP_b = OptionMenu(frBerthelot, tkP_b, *escolhasP_b, command = mudarConstante)
mnP_b.config(width = 18)

lbVU_b = Label(frBerthelot, text = "Volume")
mnV_b = OptionMenu(frBerthelot, tkV_b, *escolhasV_b, command = mudarConstante)
mnV_b.config(width = 18)

lbnU_b = Label(frBerthelot, text = "Mols")
mnn_b = OptionMenu(frBerthelot, tkn_b, *escolhasb)
mnn_b.config(width = 18)

lbTU_b = Label(frBerthelot, text = "Temperatura")
mnT_b = OptionMenu(frBerthelot, tkT_b, *escolhasT_b)
mnT_b.config(width = 18)

lbrR_b = Label(frBerthelot, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frBerthelot, text = "a")
lba_b = Label(frBerthelot, text = "m^6.Pa")
lbbb = Label(frBerthelot, text = "b")
lbb_b = Label(frBerthelot, text = "m^3.mol^-1")

btCalcularBerthelot = Button(frBerthelot, text = "Calcular", width = 18)
btCalcularBerthelot["command"] = partial(EquacaoDeBerthelotCalculo,tfP_b, tfV_b, tfn_b, lbrR_b, tfT_b, tfA_b, tfB_b, tkP_b, tkV_b)

#define the positions of the objects in the BERTHELOT
lbVazio2.grid(row = 0, column = 0)
    
lbP_b.grid(row = 1, column = 0)
tfP_b.grid(row = 1, column = 1)
lbV_b.grid(row = 2, column = 0)
tfV_b.grid(row = 2, column = 1)
lbn_b.grid(row = 3, column = 0)
tfn_b.grid(row = 3, column = 1)
lbR_b.grid(row = 4, column = 0)
lbrR_b.grid(row = 4, column = 1)
lbT_b.grid(row = 5, column = 0)
tfT_b.grid(row = 5, column = 1)
lbA_b.grid(row = 6, column = 0)
tfA_b.grid(row = 6, column = 1)
lbB_b.grid(row = 7, column = 0)
tfB_b.grid(row = 7, column = 1)

lbVazio1.grid(row = 8, column = 0)

btCalcularBerthelot.grid(row = 9, column = 1)

lbUnidades_b.grid(row = 10, column = 0)
lbPU_b.grid(row = 11, column = 0)
mnP_b.grid(row = 11, column = 1)
lbVU_b.grid(row = 12, column = 0)
mnV_b.grid(row = 12, column = 1)
lbnU_b.grid(row = 13, column = 0)
mnn_b.grid(row = 13, column = 1)
lbR_b.grid(row = 14, column = 0)
lbrR_b.grid(row = 14, column = 1)
lbba.grid(row = 15, column = 0)
lba_b.grid(row = 15, column = 1)
lbbb.grid(row = 16, column = 0)
lbb_b.grid(row = 16, column = 1)
lbTU_b.grid(row = 17, column = 0)
mnT_b.grid(row = 17, column = 1)

#define the objects of VAN DER WAALS equation
frVanDerWaals = Frame(janela)

lbVazio = Label(frVanDerWaals)
lbVazio1 = Label(frVanDerWaals)
lbVazio2 = Label(frVanDerWaals)

lbP_vdw = Label(frVanDerWaals, text = "P")
tfP_vdw = Entry(frVanDerWaals)

lbV_vdw = Label(frVanDerWaals, text = "V")
tfV_vdw = Entry(frVanDerWaals)

lbn_vdw = Label(frVanDerWaals, text = "n")
tfn_vdw = Entry(frVanDerWaals)

lbR_vdw = Label(frVanDerWaals, text = "R")

lbT_vdw = Label(frVanDerWaals, text = "T")
tfT_vdw = Entry(frVanDerWaals)

lbA_vdw = Label(frVanDerWaals, text = "a")
tfA_vdw = Entry(frVanDerWaals)

lbB_vdw = Label(frVanDerWaals, text = "b")
tfB_vdw = Entry(frVanDerWaals)

tkP_vdw = StringVar(frVanDerWaals)
escolhasP_vdw = {"Pa", "atm", "mmHg"}
tkP_vdw.set("Pa")
    
tkV_vdw = StringVar(frVanDerWaals)
escolhasV_vdw = {"m^3", "L", "mL"}
tkV_vdw.set("m^3")

tkn_vdw = StringVar(frVanDerWaals)
escolhasvdw = {"mol"}
tkn_vdw.set("mol")
    
tkT_vdw = StringVar(frVanDerWaals)
escolhasT_vdw = {"K"}
tkT_vdw.set("K")

lbUnidades_vdw = Label(frVanDerWaals, text = "Selecione as unidades da entrada")

lbPU_vdw = Label(frVanDerWaals, text = "Pressão")
mnP_vdw = OptionMenu(frVanDerWaals, tkP_vdw, *escolhasP_vdw, command = mudarConstante)
mnP_vdw.config(width = 18)

lbVU_vdw = Label(frVanDerWaals, text = "Volume")
mnV_vdw = OptionMenu(frVanDerWaals, tkV_vdw, *escolhasV_vdw, command = mudarConstante)
mnV_vdw.config(width = 18)

lbnU_vdw = Label(frVanDerWaals, text = "Mols")
mnn_vdw = OptionMenu(frVanDerWaals, tkn_vdw, *escolhasvdw)
mnn_vdw.config(width = 18)

lbTU_vdw = Label(frVanDerWaals, text = "Temperatura")
mnT_vdw = OptionMenu(frVanDerWaals, tkT_vdw, *escolhasT_vdw)
mnT_vdw.config(width = 18)

lbrR = Label(frVanDerWaals, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frVanDerWaals, text = "a")
lba = Label(frVanDerWaals, text = "m^6.Pa")
lbbb = Label(frVanDerWaals, text = "b")
lbb = Label(frVanDerWaals, text = "m^3.mol^-1")

btCalcularVanDerWaals = Button(frVanDerWaals, text = "Calcular", width = 18)
btCalcularVanDerWaals["command"] = partial(EquacaoDeVanDerWaalsCalculo,tfP_vdw, tfV_vdw, tfn_vdw, lbrR, tfT_vdw, tfA_vdw, tfB_vdw, tkP_vdw, tkV_vdw)

#define the positions of the objects in the VAN DER WAALS frame
lbVazio2.grid(row = 0, column = 0)
    
lbP_vdw.grid(row = 1, column = 0)
tfP_vdw.grid(row = 1, column = 1)
lbV_vdw.grid(row = 2, column = 0)
tfV_vdw.grid(row = 2, column = 1)
lbn_vdw.grid(row = 3, column = 0)
tfn_vdw.grid(row = 3, column = 1)
lbR_vdw.grid(row = 4, column = 0)
lbrR.grid(row = 4, column = 1)
lbT_vdw.grid(row = 5, column = 0)
tfT_vdw.grid(row = 5, column = 1)
lbA_vdw.grid(row = 6, column = 0)
tfA_vdw.grid(row = 6, column = 1)
lbB_vdw.grid(row = 7, column = 0)
tfB_vdw.grid(row = 7, column = 1)

lbVazio1.grid(row = 8, column = 0)

btCalcularVanDerWaals.grid(row = 9, column = 1)

lbUnidades_vdw.grid(row = 10, column = 0)
lbPU_vdw.grid(row = 11, column = 0)
mnP_vdw.grid(row = 11, column = 1)
lbVU_vdw.grid(row = 12, column = 0)
mnV_vdw.grid(row = 12, column = 1)
lbnU_vdw.grid(row = 13, column = 0)
mnn_vdw.grid(row = 13, column = 1)
lbR_vdw.grid(row = 14, column = 0)
lbrR.grid(row = 14, column = 1)
lbba.grid(row = 15, column = 0)
lba.grid(row = 15, column = 1)
lbbb.grid(row = 16, column = 0)
lbb.grid(row = 16, column = 1)
lbTU_vdw.grid(row = 17, column = 0)
mnT_vdw.grid(row = 17, column = 1)

#define the objects of COMPRESSIBILITY FACTOR
frFatorCompressibilidade = Frame(janela)

lbVazio = Label(frFatorCompressibilidade)
lbVazio1 = Label(frFatorCompressibilidade)
lbVazio2 = Label(frFatorCompressibilidade)

lbn_fc = Label(frFatorCompressibilidade, text = "n")
tfn_fc = Entry(frFatorCompressibilidade)

lbR_fc = Label(frFatorCompressibilidade, text = "R")

lbT_fc = Label(frFatorCompressibilidade, text = "T")
tfT_fc = Entry(frFatorCompressibilidade)

lbA_fc = Label(frFatorCompressibilidade, text = "a")
tfA_fc = Entry(frFatorCompressibilidade)

lbB_fc = Label(frFatorCompressibilidade, text = "b")
tfB_fc = Entry(frFatorCompressibilidade)

tkP_fc = StringVar(frFatorCompressibilidade)
escolhasP_fc = {"Pa", "atm", "mmHg"}
tkP_fc.set("Pa")
    
tkV_fc = StringVar(frFatorCompressibilidade)
escolhasV_fc = {"m^3", "L", "mL"}
tkV_fc.set("m^3")

tkn_fc = StringVar(frFatorCompressibilidade)
escolhasfc = {"mol"}
tkn_fc.set("mol")
    
tkT_fc = StringVar(frFatorCompressibilidade)
escolhasT_fc = {"K"}
tkT_fc.set("K")

lbUnidades_fc = Label(frFatorCompressibilidade, text = "Selecione as unidades da entrada")

lbPU_fc = Label(frFatorCompressibilidade, text = "Pressão")
mnP_fc = OptionMenu(frFatorCompressibilidade, tkP_fc, *escolhasP_fc, command = mudarConstante)
mnP_fc.config(width = 18)

lbVU_fc = Label(frFatorCompressibilidade, text = "Volume")
mnV_fc = OptionMenu(frFatorCompressibilidade, tkV_fc, *escolhasV_fc, command = mudarConstante)
mnV_fc.config(width = 18)

lbnU_fc = Label(frFatorCompressibilidade, text = "Mols")
mnn_fc = OptionMenu(frFatorCompressibilidade, tkn_fc, *escolhasfc)
mnn_fc.config(width = 18)

lbTU_fc = Label(frFatorCompressibilidade, text = "Temperatura")
mnT_fc = OptionMenu(frFatorCompressibilidade, tkT_fc, *escolhasT_fc)
mnT_fc.config(width = 18)

lbrR_fc = Label(frFatorCompressibilidade, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frFatorCompressibilidade, text = "a")
lba_fc = Label(frFatorCompressibilidade, text = "m^6.Pa")
lbbb = Label(frFatorCompressibilidade, text = "b")
lbb_fc = Label(frFatorCompressibilidade, text = "m^3.mol^-1")

btGerarGraficoFatorCompressibilidade = Button(frFatorCompressibilidade, text = "Calcular", width = 18)
btGerarGraficoFatorCompressibilidade["command"] = partial(GerarFatoresCompressibilidade, tfn_fc, lbrR_fc, tfT_fc, tfA_fc, tfB_fc, tkV_fc)

#define the positions of the objects in the COMPRESSIBILITY FACTOR frame
lbVazio2.grid(row = 0, column = 0)
    
lbn_fc.grid(row = 1, column = 0)
tfn_fc.grid(row = 1, column = 1)
lbR_fc.grid(row = 2, column = 0)
lbrR_fc.grid(row = 2, column = 1)
lbT_fc.grid(row = 3, column = 0)
tfT_fc.grid(row = 3, column = 1)
lbA_fc.grid(row = 4, column = 0)
tfA_fc.grid(row = 4, column = 1)
lbB_fc.grid(row = 5, column = 0)
tfB_fc.grid(row = 5, column = 1)

lbVazio1.grid(row = 6, column = 0)

btGerarGraficoFatorCompressibilidade.grid(row = 7, column = 1)

lbUnidades_fc.grid(row = 8, column = 0)
lbPU_fc.grid(row = 9, column = 0)
mnP_fc.grid(row = 9, column = 1)
lbVU_fc.grid(row = 10, column = 0)
mnV_fc.grid(row = 10, column = 1)
lbnU_fc.grid(row = 11, column = 0)
mnn_fc.grid(row = 11, column = 1)
lbR_fc.grid(row = 12, column = 0)
lbrR_fc.grid(row = 12, column = 1)
lbTU_fc.grid(row = 13, column = 0)
mnT_fc.grid(row = 13, column = 1)
lbba.grid(row = 14, column = 0)
lba_fc.grid(row = 14, column = 1)
lbbb.grid(row = 15, column = 0)
lbb_fc.grid(row = 15, column = 1)

#define the objects of ISOTHERMS
frIsotermas = Frame(janela)

lbVazio = Label(frIsotermas)
lbVazio1 = Label(frIsotermas)
lbVazio2 = Label(frIsotermas)

lbn_isot = Label(frIsotermas, text = "n")
tfn_isot = Entry(frIsotermas)

lbR_isot = Label(frIsotermas, text = "R")

lbT_isot = Label(frIsotermas, text = "T")
tfT_isot = Entry(frIsotermas)

lbA_isot = Label(frIsotermas, text = "a")
tfA_isot = Entry(frIsotermas)

lbB_isot = Label(frIsotermas, text = "b")
tfB_isot = Entry(frIsotermas)

tkP_isot = StringVar(frIsotermas)
escolhasP_isot = {"Pa", "atm", "mmHg"}
tkP_isot.set("Pa")
    
tkV_isot = StringVar(frIsotermas)
escolhasV_isot = {"m^3", "L", "mL"}
tkV_isot.set("m^3")

tkn_isot = StringVar(frIsotermas)
escolhasiso = {"mol"}
tkn_isot.set("mol")
    
tkT_isot = StringVar(frIsotermas)
escolhasT_isot = {"K"}
tkT_isot.set("K")

lbUnidades_isot = Label(frIsotermas, text = "Selecione as unidades da entrada")

lbPU_isot = Label(frIsotermas, text = "Pressão")
mnP_isot = OptionMenu(frIsotermas, tkP_isot, *escolhasP_isot, command = mudarConstante)
mnP_isot.config(width = 18)

lbVU_isot = Label(frIsotermas, text = "Volume")
mnV_isot = OptionMenu(frIsotermas, tkV_isot, *escolhasV_isot, command = mudarConstante)
mnV_isot.config(width = 18)

lbnU_isot = Label(frIsotermas, text = "Mols")
mnn_isot = OptionMenu(frIsotermas, tkn_isot, *escolhasiso)
mnn_isot.config(width = 18)

lbTU_isot = Label(frIsotermas, text = "Temperatura")
mnT_isot = OptionMenu(frIsotermas, tkT_isot, *escolhasT_isot)
mnT_isot.config(width = 18)

lbrR_isot = Label(frIsotermas, text = "8.31 J.K^-1.mol^-1")

lbba = Label(frIsotermas, text = "a")
lba_isot = Label(frIsotermas, text = "m^6.Pa")
lbbb = Label(frIsotermas, text = "b")
lbb_isot = Label(frIsotermas, text = "m^3.mol^-1")

btGerarGraficoIsotermas = Button(frIsotermas, text = "Gerar", width = 18)
btGerarGraficoIsotermas["command"] = partial(GerarIsotermas, tfn_isot, lbrR_isot, tfT_isot, tfA_isot, tfB_isot)

#define the positions of the objects in the ISOTHERMS frame
lbVazio2.grid(row = 0, column = 0)
    
lbn_isot.grid(row = 1, column = 0)
tfn_isot.grid(row = 1, column = 1)
lbR_isot.grid(row = 2, column = 0)
lbrR_isot.grid(row = 2, column = 1)
lbT_isot.grid(row = 3, column = 0)
tfT_isot.grid(row = 3, column = 1)
lbA_isot.grid(row = 4, column = 0)
tfA_isot.grid(row = 4, column = 1)
lbB_isot.grid(row = 5, column = 0)
tfB_isot.grid(row = 5, column = 1)

lbVazio1.grid(row = 6, column = 0)

btGerarGraficoIsotermas.grid(row = 7, column = 1)

lbUnidades_isot.grid(row = 8, column = 0)
lbPU_isot.grid(row = 9, column = 0)
mnP_isot.grid(row = 9, column = 1)
lbVU_isot.grid(row = 10, column = 0)
mnV_isot.grid(row = 10, column = 1)
lbnU_isot.grid(row = 11, column = 0)
mnn_isot.grid(row = 11, column = 1)
lbR_isot.grid(row = 12, column = 0)
lbrR_isot.grid(row = 12, column = 1)
lbTU_isot.grid(row = 13, column = 0)
mnT_isot.grid(row = 13, column = 1)
lbba.grid(row = 14, column = 0)
lba_isot.grid(row = 14, column = 1)
lbbb.grid(row = 15, column = 0)
lbb_isot.grid(row = 15, column = 1)

#Toolbar and Label on the bot of the frame with the results and warnings
frToolBarBottom = Frame(janela, bg = "green")
lbAviso = Label(frToolBarBottom, text = "Entrada inválida", bg = "red")
lbAviso.grid(row = 0, column = 0)

#Toolbar with the equations and transformations
frToolBar = Frame(janela, bg = "gray")

#array to change save the previous screens
global vetor
vetor = [""]

btFatoresCompressibilidade = Button(frToolBar, text = "Fatores de Compressibilidade", width = 25)
btFatoresCompressibilidade["command"] = partial(escolhaTela, "FatoresCompressibilidade")
btIsotermas = Button(frToolBar, text = "Isotermas", width = 25)
btIsotermas["command"] = partial(escolhaTela, "Isotermas")

tkEquacoesDeEstado = StringVar(frToolBar)
escolhasEquacoesDeEstado = {"van der Waals", "Berthelot", "Dieterici", "Redlich-Kwong"}
tkEquacoesDeEstado.set("van der Waals")

mnEquacoesDeEstado = OptionMenu(frToolBar, tkEquacoesDeEstado, *escolhasEquacoesDeEstado, command = escolhaTela) #OLHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
mnEquacoesDeEstado["width"] = 25

#put the buttons on the toolbar
mnEquacoesDeEstado.grid(row = 0, column = 0, padx = 5)
btFatoresCompressibilidade.grid(row = 0, column = 1, padx = 5)
btIsotermas.grid(row = 0, column = 2, padx = 5)

#put the toolbar on the mainframe
frToolBar.grid(row = 0, column = 0, sticky = W+E)

janela.mainloop()