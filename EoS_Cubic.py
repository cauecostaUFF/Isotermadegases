import numpy as np
import matplotlib.pyplot as plt
# import streamlit as st

class Component:
    def __init__(self,name =' ',molarmass=0,Tc=0,Pc=0,Vc=0,Zc=0,omega=0,Ttriple=0,Ptriple=0):
        # Definindo os campos
        self.Name = name
        self.MolarMass = molarmass
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Zc = Zc
        self.omega = omega
        # Definindo dados experimentais de Psat
        self.dados = np.zeros(1)
        self.Ttriple = Ttriple
        self.Ptriple = Ptriple

    def Rackett(self,T):
        vol = self.Vc*self.Zc**((1-T/self.Tc)**(2/7))
        return vol


class Cubic_eos:
    def __init__(self,nome):
        self.Name = nome
        if nome=='vdW':
            sig = 0
            e = 0
            om = 1 / 8
            psi = 27 / 64
            zc = 3 / 8
            alp = lambda Tr, omega: 1
        elif nome == "RK":
            sig = 1
            e = 0
            om = .08664
            psi = .42748
            zc = 1 / 3
            alp = lambda Tr, omega: 1/np.sqrt(Tr)
        elif nome == "SRK":
            sig = 1
            e = 0
            om = .08664
            psi = .42748
            zc = 1 / 3
            alp = lambda Tr,omega: (1+(.480+1.574*omega-.176*omega**2)*(1-Tr**.5))**2
        elif nome == "PR":
            sig=1 + np.sqrt(2)
            e= 1 - np.sqrt(2)
            om=.07780
            psi=.45724
            zc=.30740
            alp = lambda Tr, omega: (1 + (.37464+1.54226*omega-.26992*omega**2) * (1 - Tr ** .5)) ** 2
        self.sigma = sig;
        self.ep = e;
        self.Omega = om;
        self.Psi = psi;
        self.Zc = zc;
        self.alpha = alp
        self.R = 83.14

    def CalcP(eos,comp,T,V):
        a = eos.Psi * eos.alpha(T / comp.Tc, comp.omega)* eos.R**2. * comp.Tc**2. / comp.Pc
        b = eos.Omega * eos.R * comp.Tc / comp.Pc
        P = eos.R * T / (V - b) - a / (V + eos.ep * b) / (V + eos.sigma * b)
        return P

    def CalcPhi(eos,comp,T,V):
            a = eos.Psi*eos.alpha(T/comp.Tc,comp.omega)*eos.R**2.*comp.Tc**2./comp.Pc
            b = eos.Omega*eos.R*comp.Tc/comp.Pc
            P = eos.CalcP(comp,T,V)
            beta = b*P/eos.R/T
            Z = P*V/eos.R/T
            if eos.Name=='vdW':
                I = beta/Z
            else:
                I = 1/(eos.sigma-eos.ep) * np.log((Z+eos.sigma*beta)/(Z+eos.ep*beta))

            q_ =  eos.Psi*eos.alpha(T/comp.Tc,comp.omega)/eos.Omega/(T/comp.Tc)
            lnphi = (Z-1) - np.log(Z-beta) - q_*I
            return lnphi,P

    def EoSRoots(eos,comp,T,P):
        alp = eos.alpha(T/comp.Tc,comp.omega)
        A= eos.Psi*alp*eos.R**2*comp.Tc**2/comp.Pc
        B= eos.Omega*eos.R*comp.Tc/comp.Pc
        CA = A*P/eos.R**2/T**2
        CB = P*B/eos.R/T
        if eos.Name=="vdW":
            Vec = [1, -(1 + CB), CA, -CA*CB]
        elif eos.Name == "PR":
            Vec = [1, -(1 - CB),(CA -3*CB**2 - 2*CB), -(CA*CB - CB**2 - CB**3)]
        else:
            Vec = [1,-1, (CA-CB -CB**2), -CA*CB]
        Z = np.roots(Vec)
        V = Z * eos.R * T / P
        return V

    def ELVPure(eos,comp,T):

        P0 = comp.Pc * np.exp(5.4 * (comp.omega + 1) * (1 - comp.Tc / T))
        Vol = eos.EoSRoots(comp, T, P0)

        Vliq = np.min(Vol)
        Vvap = np.max(Vol)

        phi_liq,P = eos.CalcPhi(comp, T, Vliq)
        phi_vap,P = eos.CalcPhi(comp, T, Vvap)

        phi = np.array([phi_liq,phi_vap])
        cont = 0
        P1 = P0
        P0 = (P1 * np.exp(phi[0]) / np.exp(phi[1]))

        while abs(np.exp(phi[0]) / np.exp(phi[1]) - 1) > 1e-4:
            Vol = eos.EoSRoots(comp, T, P0)

            Vliq = min(Vol)
            Vvap = max(Vol)

            phi_liq,P = eos.CalcPhi(comp, T, Vliq)
            phi_vap,P = eos.CalcPhi(comp, T, Vvap)
            phi = np.array([phi_liq, phi_vap])
            cont += 1
            P1 = P0
            P0 = np.real(P1 * np.exp(phi[0]) / np.exp(phi[1]))

        return P0,[Vliq,Vvap]

    def PlotGraf_PV(eos,comp,T1=0):
        T = np.linspace(comp.Ttriple,comp.Tc,101)
        P = np.zeros(len(T))
        V = np.zeros([len(T),2])
        for i in range(len(T)):
            [P[i], V[i]] = eos.ELVPure(comp,T[i])

        if T1 == 0:
            T1 = comp.Tc

        Vr = np.logspace(np.log10(np.min(V)),np.log10(np.max(V)),1001)
        Pr = np.zeros(len(Vr))
        Pc,Vc = eos.ELVPure(comp, comp.Tc)
        for i in range(len(Vr)):
            if T1 <= comp.Tc:
                Peq, Veq = eos.ELVPure(comp, T1)
                if Vr[i]<np.min(Veq) or Vr[i]>np.max(Veq):
                    Pr[i] = eos.CalcP(comp,T1,Vr[i])
                else:
                    Pr[i] = Peq
            else:
                Pr[i] = eos.CalcP(comp, T1, Vr[i])
        p = plt.loglog(V,P,'k',Vc[0],Pc,'or',Vr,Pr,markersize=4)
        plt.xlabel(r'Volume molar (cm$^3$/mol)')
        plt.ylabel('Pressão (bar)')
        plt.legend([p[0],p[3],p[2]],['Curva de equilíbrio','Isoterma à '+str(np.round(T1))+'K','Ponto crítico'])
        plt.show()

# eos = Cubic_eos('PR')

# comp = Component('Water',1.801530e+01,6.471300e+02,2.205500e+02,5.594780e+01,2.290000e-01,3.448610e-01,273.16,0.006)

# lnphi,P = eos.CalcPhi(comp,300,30)
# print(lnphi,P)

# V = eos.EoSRoots(comp,300,1)
# print(V)
# P = eos.ELVPure(comp,373.15)
# print(P)
# eos.PlotGraf_PV(comp,680)

