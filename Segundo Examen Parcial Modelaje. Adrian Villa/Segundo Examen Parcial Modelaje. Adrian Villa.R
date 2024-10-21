
## ADRIÁN VILLA ##

## SEGUNDO EXAMEN PARCIAL MODELAJE DE ENFERMEDADES INFECCIOSAS ##
# 20 / 10 / 2024 #

# PROBLEMA 1. SELECCIÓN.

# PROBLEMA 2. ESTADIO ASINTOMÁTICO

# 1. Dibuja el esquema por compartimentos que representa este conjunto de ecuaciones


# 2. Explica el significado de cada ecuación: es decir, explica el significado de cada término y cada parámetro. 
#    ¿Existen parámetros que están acotados, si es así indica el intervalo en donde pueden variar?

# dS = upsilon - beta * S * (I + q * A) - mu * S
# dE = beta * S * (I + q * A) - (eta + mu) * E
# dI = p * eta * E - (alpha + mu) * I
# dA = (1 - p) * eta * E - (gamma + mu) * A
# dR = alpha * I + gamma * A - mu * R

# S = Susceptibles, E = Expuestos, I = Infectados (sintomáticos), A = Asintomáticos (infectivos), R = Recuperados
# beta es la tasa de infección, beta * q es la tasa de infección "ajustada" de los asintomáticos infectivos
# eta es la tasa de progreso del patógeno en un individio expuesto
# p * eta * E es la proporción de Expuestos que se vuelven infectados
# gamma es la tasa de recuperación de los Asintomáticos
# alpha es la tasa de recuperación de los Infectados
# mu es la tasa de muerte de cada grupo
# upsilon es la tasa de nacimiento

# p y q están acotados, sus valores deben ser entre 0 y 1
# beta, eta, gamma, alpha, mu, upsilon son tasas no negativas, deben ser iguales o mayores a 0

# 3. ¿Bajo qué condiciones la población se conservaría?
# (N(t) = S(t) + E(t) + I(t) + A(t) + R(t)) / dt
# Sustituimos cada ecuación
# N(t) = upsilom - mu * N
# upsilon = N * mu -> En resumen, la tasa de natalidad tiene que ser igual a la suma de la tasa de mortalidad de todos los grupos

# 4. Encuentra, si existe, el punto de equilibrio free-disease
# Free-disease quiere decir que toda la población es susceptible por lo que tenemos un p.e (?, 0, 0, 0, 0)
# Como toda la población está contenida en S, sustiimos los valores y queda dS/dt = upsilon - mu * S
# upsilon - mu * S = 0, despejamos y S = upsilon / mu
# P.E. Free-Disease = (Upsilon / mu, 0, 0, 0, 0)

# 5. ¿Qué tipo de enfermedad puede estar describiendo? Justifica tu respuesta.
# Lo clásico sería pensar en COVID-19 (o la mayoría de sus variantes), lo cual se ajusta bastante bien al modelo
# Otra enfermedad podría ser la poliomielitis, cerca del 70% de los casos son asintomáticos, pero estos casos juegan un
#  papel importante a la hora de contrastar la dinámica infecciosa.

# 6. Selecciona un conjunto de parámetros adecuados y resuelve numéricamente el sistema de ecuaciones diferenciales. 
#    Asegurate que tu solución alcance un punto de equilibrio. Discute tu resultado.
# Instalar librerías necesarias
install.packages("deSolve")
library (deSolve)
install.packages("dplyr")
library (dplyr)
install.packages("nleqslv")
library(nleqslv)

# Definir los parámetros del sistema
beta <- 0.8  # Tasa de transmisión
q <- 0.9     # Proporción de la transmisibilidad de asintomáticos
eta <- 0.4   # Tasa de progresión de expuestos a infecciosos
p <- 0.6     # Proporción de los expuestos que se convierten en infectados sintomáticos
alpha <- 0.1 # Tasa de recuperación de infectados sintomáticos
gamma <- 0.1 # Tasa de recuperación de asintomáticos
mu <- 0.05   # Tasa de mortalidad natural
nu <- 0.08   # Tasa de natalidad

# Función que representa el sistema en equilibrio (es decir, con las derivadas igual a 0)
seiar_equilibrium <- function(x) {
  S <- x[1]  # Susceptibles
  E <- x[2]  # Expuestos
  I <- x[3]  # Infectados sintomáticos
  A <- x[4]  # Asintomáticos
  R <- x[5]  # Recuperados
  
  # Ecuaciones en equilibrio (dS/dt = dE/dt = dI/dt = dA/dt = dR/dt = 0)
  eq1 <- nu - beta * S * (I + q * A) - mu * S              # dS/dt = 0
  eq2 <- beta * S * (I + q * A) - (eta + mu) * E            # dE/dt = 0
  eq3 <- p * eta * E - (alpha + mu) * I                     # dI/dt = 0
  eq4 <- (1 - p) * eta * E - (gamma + mu) * A               # dA/dt = 0
  eq5 <- alpha * I + gamma * A - mu * R                     # dR/dt = 0
  
  return(c(eq1, eq2, eq3, eq4, eq5))
}

# Valores iniciales para S, E, I, A, R
initial_values <- c(S = 980, E = 10, I = 5, A = 5, R = 0)

# Resolver el sistema no lineal usando nleqslv
result <- nleqslv(initial_values, seiar_equilibrium)

# Ver los resultados del punto de equilibrio
result$x  # Valores de S, E, I, A, R en el punto de equilibrio

# Si modelamos el Free-disease (que solo haya población en S), el P.E solo nos dará valores para S
# Si cambiamos los parámetros (añadiendo E, I, A) obtenemos valores numéricos para cada sector, este valor viene en %
# Lo que yo obtuve con estos parámetros:    S         E         I         A         R 
#                                       0.2197266 0.1533637 0.2453819 0.1635880 0.8179398 
