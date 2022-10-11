import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import EoS_Cubic as EoS
# from matplotlib.widgets import Button
# from matplotlib.widgets import TextBox

st.session_state.visibility = "collapsed"
st.session_state.disabled = True
st.set_page_config(
    page_title="Isoterma de gases",
    # page_icon="🧊",
    layout="centered",
    initial_sidebar_state="expanded"
)

# Definindo modelos
modelos = {0: r"van der Waals", 1: r"Redlich-Kwong", 2: r"Soave-Redlich-Kwong", 3: r"Peng-Robinson"}
def format_func(option):
    return modelos[option]
nome = ['vdW','RK','SRK','PR']
# Definindo substâncias e suas propriedades
subs = {0:"Água",
        1:"Metano",
        2:"Gás Carbônico",
        3:"Ciclohexano",
        4:"Benzeno",
        5:"Metanol",
        6:"Etanol"}
Tcval = np.array([6.471300e+02,1.905640e+02,3.0421e2,5.538000e+02,5.620500e+02,5.125000e+02,5.140000e+02])
Pcval = np.array([2.205500e+02,4.599000e+01,7.383e1,4.080000e+01,4.895000e+01,8.084000e+01,6.137000e+01])
omegaval = np.array([3.448610e-01,1.154780e-02,2.236210e-01,2.080540e-01,2.103000e-01,5.658310e-01,6.435580e-01])
Ttriple = np.array([273.15,9.069000e+01,216.55,2.796900e+02,2.786800e+02,1.754700e+02,1.590500e+02])
def format_subs(option):
    return subs[option]

# Informações
extra, main_col1,main_col2,main_col3 = st.columns([0.5,10,10,10])

with main_col2:
    st.write(r"Selecione a substância:")
    option_subs = st.selectbox(
        "Água",
        options=list(subs.keys()), format_func=format_subs,
        label_visibility=st.session_state.visibility,
        disabled=False,
    )

with main_col1:
    # Selecionando modelo
    st.write(r"Selecione a equação de estado:")
    # st.write(r"Modelo de Margules: G$^{E}$ = n$\beta$RTx$_A$x$_B$")
    option_eq = st.selectbox(
        "van der Waals",
        options=list(modelos.keys()),format_func=format_func,
        label_visibility=st.session_state.visibility,
        disabled=False,
    )
with main_col3:
    st.write(r"Defina a Temperatura (em K):")
    T1 = st.number_input("", min_value=Ttriple[option_subs], max_value=1000.0, value=Tcval[option_subs]*.75, step=2.5,
                         label_visibility=st.session_state.visibility)

extra, lin1_col1,lin1_col2,lin1_col3 = st.columns([0.5,10,10,10])
with lin1_col3:
    st.write(r"Fator acêntrico: ")
    st.write(str(omegaval[option_subs]))
    # omega = st.number_input("", min_value=0.0,max_value=1.0,value = 0.229,step=0.001,format="%.3f", label_visibility=st.session_state.visibility)

with lin1_col2:
    st.write(r"Pressão crítica: ")
    st.write(str(Pcval[option_subs]) + " bar")
    # Pc = st.number_input("", min_value=0.0,max_value=1000.0,value = Pcval[option_subs] ,step=1.0, label_visibility=st.session_state.visibility)
with lin1_col1:
    st.write(r"Temperatura crítica: ")
    st.write(str(Tcval[option_subs])+" K")
    # Tc = st.number_input("", min_value=0.0,max_value=1000.0,value =Tcval[option_subs] ,step=2.5, label_visibility=st.session_state.visibility)

eos = EoS.Cubic_eos(nome[option_eq])
# comp = EoS.Component(format_func(option_eq),Tc,Pc,Vc = 100,Zc = eos.Zc ,omega=omega,Ttriple=273)
comp = EoS.Component(format_subs(option_subs),0,Tcval[option_subs],Pcval[option_subs],0,0,omegaval[option_subs],Ttriple[option_subs],0)
# st.write(Tc,Pc,omega,eos.Name)
eos.PlotGraf_PV(comp,T1)
# P = eos.ELVPure(comp,373.15)
# st.write(P)
plt.show()
st.pyplot(plt)
col_end1,col_end2,col_end3 = st.columns(3)
with col_end3:
    st.write("Desenvolvido por: Cauê Costa\nE-mail: cauecosta@id.uff.br")
with col_end2:
    st.image("UFF.png", width=150)
with col_end1:
    st.image("molmod.png", width=150)