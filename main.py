# =========================================================
# ORGANICCHEMAI PRO MAX
# VERSION UNIVERSITARIA DEFINITIVA
# =========================================================

import streamlit as st
import pandas as pd
import re

# =========================================================
# CONFIGURACIÓN
# =========================================================

st.set_page_config(
    page_title="OrganicChemAI PRO MAX",
    page_icon="🧪",
    layout="wide"
)

# =========================================================
# ESTILOS
# =========================================================

st.markdown("""
<style>

html, body, [class*="css"]{
    background-color:#0b1120;
    color:white;
    font-family:'Segoe UI';
}

.main{
    background-color:#0b1120;
}

.title-box{
    background:linear-gradient(145deg,#111827,#1f2937);
    padding:35px;
    border-radius:25px;
    margin-bottom:25px;
    border:1px solid rgba(255,255,255,0.08);
}

.card{
    background:#111827;
    padding:22px;
    border-radius:20px;
    margin-bottom:20px;
    border-left:5px solid #00ffd5;
    box-shadow:0px 0px 15px rgba(0,255,213,0.08);
}

.info-box{
    background:#1e293b;
    padding:18px;
    border-radius:18px;
    margin-bottom:15px;
}

.highlight{
    color:#00ffd5;
    font-weight:bold;
}

.small-title{
    color:#00ffd5;
    font-size:20px;
    font-weight:bold;
}

</style>
""", unsafe_allow_html=True)

# =========================================================
# BASE DE DATOS IA
# =========================================================

base_datos = {

    "etanol":{
        "tipo":"Alcohol primario",
        "grupo":"Hidroxilo (-OH)",
        "polaridad":"Alta",
        "reactividad":"Oxidación y deshidratación",
        "uso":"Combustibles, bebidas alcohólicas y laboratorio.",
        "dato":"Forma puentes de hidrógeno.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/ethanol/PNG"
    },

    "metanol":{
        "tipo":"Alcohol",
        "grupo":"Hidroxilo (-OH)",
        "polaridad":"Alta",
        "reactividad":"Muy inflamable",
        "uso":"Disolvente y producción industrial.",
        "dato":"Altamente tóxico para humanos.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/methanol/PNG"
    },

    "acetona":{
        "tipo":"Cetona",
        "grupo":"Carbonilo",
        "polaridad":"Media",
        "reactividad":"Adición nucleofílica",
        "uso":"Solvente industrial y cosmético.",
        "dato":"Evapora rápidamente.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetone/PNG"
    },

    "benceno":{
        "tipo":"Compuesto aromático",
        "grupo":"Anillo aromático",
        "polaridad":"Baja",
        "reactividad":"Sustitución electrofílica",
        "uso":"Petroquímica y síntesis.",
        "dato":"Presenta resonancia electrónica.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/benzene/PNG"
    },

    "fenol":{
        "tipo":"Fenol aromático",
        "grupo":"OH aromático",
        "polaridad":"Media",
        "reactividad":"Oxidación y sustitución",
        "uso":"Medicamentos y resinas.",
        "dato":"Más ácido que alcoholes comunes.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/phenol/PNG"
    },

    "acido acetico":{
        "tipo":"Ácido carboxílico",
        "grupo":"COOH",
        "polaridad":"Alta",
        "reactividad":"Ácido débil",
        "uso":"Vinagre y síntesis.",
        "dato":"Responsable del olor del vinagre.",
        "imagen":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetic%20acid/PNG"
    }

}

# =========================================================
# ANALIZADOR
# =========================================================

def analizar_formula(formula):

    try:

        formula = formula.upper()

        C = re.search(r'C(\d*)', formula)
        H = re.search(r'H(\d*)', formula)
        O = re.search(r'O(\d*)', formula)
        N = re.search(r'N(\d*)', formula)

        c = int(C.group(1)) if C and C.group(1) else 1 if C else 0
        h = int(H.group(1)) if H and H.group(1) else 1 if H else 0
        o = int(O.group(1)) if O and O.group(1) else 1 if O else 0
        n = int(N.group(1)) if N and N.group(1) else 1 if N else 0

        insaturacion = ((2*c)+2-h+n)/2

        interpretacion = ""

        if insaturacion == 0:
            interpretacion = "La molécula parece saturada, característica típica de alcanos o cadenas simples."

        elif insaturacion >= 1 and insaturacion < 4:
            interpretacion = "La estructura podría contener dobles enlaces, ciclos o insaturaciones moderadas."

        elif insaturacion >= 4:
            interpretacion = "Existe una alta probabilidad de aromaticidad o múltiples insaturaciones."

        if o >= 1:
            interpretacion += " Además, contiene oxígeno, por lo que podría pertenecer a alcoholes, cetonas o ácidos."

        if n >= 1:
            interpretacion += " También presenta nitrógeno, indicando posible presencia de aminas o amidas."

        return {
            "C":c,
            "H":h,
            "O":o,
            "N":n,
            "I":round(insaturacion,2),
            "interpretacion":interpretacion
        }

    except:
        return None

# =========================================================
# GENERADOR IUPAC
# =========================================================

def generar_iupac(formula):

    try:

        formula = formula.upper()

        C = re.search(r'C(\d*)', formula)
        H = re.search(r'H(\d*)', formula)
        O = re.search(r'O(\d*)', formula)
        N = re.search(r'N(\d*)', formula)

        c = int(C.group(1)) if C and C.group(1) else 1
        h = int(H.group(1)) if H and H.group(1) else 1

        prefijos = {
            1:"met",
            2:"et",
            3:"prop",
            4:"but",
            5:"pent",
            6:"hex",
            7:"hept",
            8:"oct",
            9:"non",
            10:"dec",
            11:"undec",
            12:"dodec"
        }

        base = prefijos.get(c,"cadena compleja")

        if O and h == (2*c)+2:
            return base + "anol"

        elif O:
            return base + "anona"

        elif N:
            return base + "anamina"

        elif h == (2*c)+2:
            return base + "ano"

        elif h == (2*c):
            return base + "eno"

        elif h == (2*c)-2:
            return base + "ino"

        else:
            return "Compuesto orgánico complejo"

    except:
        return "No reconocido"

# =========================================================
# MENÚ
# =========================================================

menu = st.sidebar.selectbox(
    "📚 Navegación",
    [
        "Inicio",
        "Analizador Molecular",
        "Consulta IA",
        "Comparador Molecular",
        "Isomería",
        "Reacciones Orgánicas",
        "Mapa Funcional",
        "Generador IUPAC",
        "Quiz Inteligente"
    ]
)

# =========================================================
# INICIO
# =========================================================

if menu == "Inicio":

    st.title("🧪 OrganicChemAI PRO MAX")

    st.markdown("""
<div class="title-box">

# Plataforma Inteligente de Química Orgánica

Aplicación diseñada para:
- análisis molecular,
- simulación conceptual,
- interpretación química,
- nomenclatura IUPAC,
- isomería,
- evaluación interactiva.

Inspirada en:
McMurry, Chang, Morrison & Boyd.

</div>
""", unsafe_allow_html=True)

    col1,col2 = st.columns(2)

    with col1:

        st.markdown("""
<div class="card">

## 🎯 Objetivo

Facilitar el aprendizaje universitario
de química orgánica mediante herramientas
interactivas y modernas.

</div>
""", unsafe_allow_html=True)

    with col2:

        st.markdown("""
<div class="card">

## 👨‍🔬 Equipo del Proyecto

- Ana Maria Florez Barbosa
- Ines Maria Rincon
- Davis Santiago Celis

</div>
""", unsafe_allow_html=True)

# =========================================================
# ANALIZADOR MOLECULAR
# =========================================================

elif menu == "Analizador Molecular":

    st.header("🔬 Analizador Molecular")

    st.markdown("""
<div class="card">

<div class="small-title">📘 ¿Qué hace esta sección?</div>

Analiza:
- saturación,
- heteroátomos,
- insaturaciones,
- aromaticidad,
- complejidad molecular.

### ✅ Ejemplos

C2H6  
C6H6  
C2H5OH  
C6H12O6  

</div>
""", unsafe_allow_html=True)

    formula = st.text_input("Introduce una fórmula molecular")

    if formula:

        datos = analizar_formula(formula)

        if datos:

            c1,c2,c3,c4,c5 = st.columns(5)

            c1.metric("Carbonos", datos["C"])
            c2.metric("Hidrógenos", datos["H"])
            c3.metric("Oxígenos", datos["O"])
            c4.metric("Nitrógenos", datos["N"])
            c5.metric("Insaturación", datos["I"])

            st.markdown(f"""
<div class="card">

## ⚛️ Interpretación Molecular

{datos['interpretacion']}

</div>
""", unsafe_allow_html=True)

# =========================================================
# CONSULTA IA
# =========================================================

elif menu == "Consulta IA":

    st.header("🧠 Consulta IA")

    st.markdown("""
<div class="card">

<div class="small-title">📘 ¿Qué hace esta IA?</div>

Analiza:
- grupo funcional,
- polaridad,
- reactividad,
- aplicaciones,
- propiedades químicas.

### ✅ Ejemplos

etanol  
metanol  
acetona  
benceno  
fenol  
acido acetico  

</div>
""", unsafe_allow_html=True)

    consulta = st.text_input("Introduce un compuesto")

    if consulta:

        consulta = consulta.lower()

        encontrado = False

        for compuesto in base_datos:

            if consulta in compuesto:

                data = base_datos[compuesto]

                col1,col2 = st.columns([1,2])

                with col1:
                    st.image(data["imagen"], width=250)

                with col2:

                    st.markdown(f"""
<div class="card">

# 🔬 {compuesto.title()}

### ⚛️ Tipo
{data['tipo']}

### 🧪 Grupo funcional
{data['grupo']}

### 💧 Polaridad
{data['polaridad']}

### 🔥 Reactividad
{data['reactividad']}

### 🌎 Aplicaciones
{data['uso']}

### 💡 Dato interesante
{data['dato']}

</div>
""", unsafe_allow_html=True)

                encontrado = True
                break

        if not encontrado:

            st.warning("⚠️ El compuesto aún no está registrado completamente en la base de datos.")

# =========================================================
# COMPARADOR
# =========================================================

elif menu == "Comparador Molecular":

    st.header("⚖️ Comparador Molecular")

    st.markdown("""
<div class="card">

<div class="small-title">📘 ¿Qué hace esta sección?</div>

Permite comparar:
- polaridad,
- estructura,
- grupos funcionales,
- comportamiento químico.

### ✅ Ejemplo

etanol VS acetona

</div>
""", unsafe_allow_html=True)

    a = st.text_input("Primer compuesto")
    b = st.text_input("Segundo compuesto")

    if a and b:

        st.markdown(f"""
<div class="card">

# 🔬 Comparación Molecular

### {a.title()} VS {b.title()}

- Diferencias estructurales importantes.
- Variaciones de polaridad y reactividad.
- Posibles diferencias en punto de ebullición y solubilidad.

</div>
""", unsafe_allow_html=True)

# =========================================================
# ISOMERÍA
# =========================================================

elif menu == "Isomería":

    st.header("🧬 Isomería Inteligente")

    st.markdown("""
<div class="card">

<div class="small-title">📘 ¿Qué analiza esta sección?</div>

Detecta posibles:
- isómeros de cadena,
- posición,
- función,
- geométricos.

### ✅ Ejemplos

C4H10  
C2H6O  
C4H8  
C5H12  

</div>
""", unsafe_allow_html=True)

    iso = st.text_input("Introduce una fórmula molecular")

    if iso:

        resultado = ""

        if iso.upper() == "C4H10":
            resultado = "Presenta isomería de cadena debido a diferentes posibles distribuciones del esqueleto carbonado."

        elif iso.upper() == "C2H6O":
            resultado = "Presenta isomería funcional: puede formar alcoholes o éteres."

        elif iso.upper() == "C4H8":
            resultado = "Puede presentar isomería geométrica y de posición."

        elif iso.upper() == "C5H12":
            resultado = "Presenta múltiples isómeros de cadena por reorganización de carbonos."

        else:
            resultado = "La fórmula podría presentar diferentes configuraciones estructurales dependiendo de la distribución molecular."

        st.markdown(f"""
<div class="card">

## 🔄 Resultado del análisis

{resultado}

</div>
""", unsafe_allow_html=True)

# =========================================================
# REACCIONES ORGÁNICAS
# =========================================================

elif menu == "Reacciones Orgánicas":

    st.header("⚗️ Reacciones Orgánicas")

    reaccion = st.selectbox(
        "Selecciona una reacción",
        [
            "Adición",
            "Sustitución",
            "Eliminación",
            "Oxidación",
            "Reducción",
            "Esterificación",
            "Hidrogenación"
        ]
    )

    if reaccion == "Adición":

        st.markdown("""
<div class="card">

# ➕ Adición

### 🧪 Ejemplo

C2H4 + H2 → C2H6

### 📚 Explicación

Los alquenos rompen enlaces π
para incorporar nuevos átomos.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Sustitución":

        st.markdown("""
<div class="card">

# 🔄 Sustitución

### 🧪 Ejemplo

CH4 + Cl2 → CH3Cl + HCl

### 📚 Explicación

Un átomo es reemplazado
por otro grupo funcional.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Eliminación":

        st.markdown("""
<div class="card">

# ⬆️ Eliminación

### 🧪 Ejemplo

C2H5OH → C2H4 + H2O

### 📚 Explicación

La molécula pierde grupos
para formar dobles enlaces.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Oxidación":

        st.markdown("""
<div class="card">

# 🔥 Oxidación

### 🧪 Ejemplo

Alcohol → Aldehído → Ácido

### 📚 Explicación

Incrementan enlaces con oxígeno.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Reducción":

        st.markdown("""
<div class="card">

# 🔬 Reducción

### 🧪 Ejemplo

Aldehído → Alcohol

### 📚 Explicación

Aumentan hidrógenos
o disminuye oxígeno.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Esterificación":

        st.markdown("""
<div class="card">

# 🧴 Esterificación

### 🧪 Ejemplo

Ácido + Alcohol → Éster + Agua

### 📚 Explicación

Forma compuestos aromáticos
y derivados industriales.

</div>
""", unsafe_allow_html=True)

    elif reaccion == "Hidrogenación":

        st.markdown("""
<div class="card">

# 💨 Hidrogenación

### 🧪 Ejemplo

Alqueno + H2 → Alcano

### 📚 Explicación

Se adiciona hidrógeno
a dobles enlaces.

</div>
""", unsafe_allow_html=True)

# =========================================================
# MAPA FUNCIONAL
# =========================================================

elif menu == "Mapa Funcional":

    st.header("🗺️ Mapa Funcional")

    grupos = []

    lista = [
        ("Alcohol","-OH","Etanol"),
        ("Cetona","C=O","Acetona"),
        ("Aldehído","-CHO","Propanal"),
        ("Ácido","-COOH","Ácido acético"),
        ("Éster","-COOR","Etanoato"),
        ("Amina","-NH2","Etilamina"),
        ("Amida","-CONH2","Acetamida"),
        ("Éter","R-O-R","Dietil éter"),
        ("Alqueno","C=C","Eteno"),
        ("Alquino","C≡C","Etino")
    ]

    for i in range(5):
        for g in lista:
            grupos.append({
                "Grupo":g[0],
                "Fórmula":g[1],
                "Ejemplo":g[2]
            })

    df = pd.DataFrame(grupos)

    st.dataframe(df)

# =========================================================
# IUPAC
# =========================================================

elif menu == "Generador IUPAC":

    st.header("🧾 Generador IUPAC")

    st.markdown("""
<div class="card">

<div class="small-title">📘 ¿Qué hace esta sección?</div>

Genera nombres aproximados
según la fórmula molecular.

### ✅ Ejemplos

C4H10  
C3H6  
C2H6O  
C5H10  

</div>
""", unsafe_allow_html=True)

    formula = st.text_input("Introduce una fórmula molecular")

    if formula:

        nombre = generar_iupac(formula)

        st.success(f"✅ Nombre generado: {nombre}")

# =========================================================
# QUIZ
# =========================================================

elif menu == "Quiz Inteligente":

    st.header("🎓 Quiz Inteligente")

    preguntas = [

        {
            "pregunta":"¿Qué grupo funcional posee el etanol?",
            "opciones":["OH","COOH","NH2"],
            "correcta":"OH",
            "tema":"Alcoholes"
        },

        {
            "pregunta":"¿Qué compuesto es aromático?",
            "opciones":["Benceno","Metano","Etanol"],
            "correcta":"Benceno",
            "tema":"Aromaticidad"
        },

        {
            "pregunta":"¿Qué reacción forma alquenos?",
            "opciones":["Eliminación","Oxidación","Reducción"],
            "correcta":"Eliminación",
            "tema":"Reacciones Orgánicas"
        },

        {
            "pregunta":"¿Qué grupo funcional poseen los ácidos carboxílicos?",
            "opciones":["COOH","OH","NH2"],
            "correcta":"COOH",
            "tema":"Ácidos carboxílicos"
        },

        {
            "pregunta":"¿Qué compuesto contiene triple enlace?",
            "opciones":["Alquino","Alcohol","Cetona"],
            "correcta":"Alquino",
            "tema":"Alquinos"
        }

    ]

    if "indice" not in st.session_state:
        st.session_state.indice = 0

    if "puntaje" not in st.session_state:
        st.session_state.puntaje = 0

    if "errores" not in st.session_state:
        st.session_state.errores = []

    if st.session_state.indice < len(preguntas):

        p = preguntas[st.session_state.indice]

        st.markdown(f"""
<div class="card">

# ❓ Pregunta {st.session_state.indice+1}/5

{p['pregunta']}

</div>
""", unsafe_allow_html=True)

        respuesta = st.radio(
            "Selecciona una opción",
            p["opciones"],
            key=st.session_state.indice
        )

        if st.button("Responder"):

            if respuesta == p["correcta"]:

                st.success("✅ Correcto")
                st.session_state.puntaje += 1

            else:

                st.error("❌ Incorrecto")
                st.session_state.errores.append(p["tema"])

            st.session_state.indice += 1
            st.rerun()

    else:

        st.markdown("""
<div class="card">

# 🏆 RESULTADO FINAL

</div>
""", unsafe_allow_html=True)

        st.write(f"## Puntaje: {st.session_state.puntaje}/5")

        if st.session_state.puntaje == 5:

            st.balloons()
            st.success("🎉 Excelente dominio de química orgánica.")

        elif st.session_state.puntaje >= 3:

            st.success("✅ Buen desempeño general.")

        else:

            st.error("⚠️ Debes reforzar varios conceptos fundamentales.")

        if len(st.session_state.errores) > 0:

            st.warning(
                "📘 Temas recomendados para estudiar: "
                + ", ".join(set(st.session_state.errores))
            )

        if st.button("Reiniciar Quiz"):

            st.session_state.indice = 0
            st.session_state.puntaje = 0
            st.session_state.errores = []
            st.rerun()
            