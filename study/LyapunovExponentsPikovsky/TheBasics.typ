#import "../preamble.typ": *

#set text(
  font: "JetBrainsMono NF",
  fill: textcolors.at(0), 
  lang: "es")
#let up = sym.arrow.t
#let down = sym.arrow.b
#show emph: set text(textcolors.at(1))
#show strong: set text(textcolors.at(2))
#show: thmrules

#set page(fill: backgroundcolor)

// #show heading.where(level: 1): set text(headingcolors.at(0))
// #show heading.where(level: 2): set text(headingcolors.at(1))
// #show heading.where(level: 3): set text(headingcolors.at(2))
// #show heading.where(level: 4): set text(headingcolors.at(3))


#set math.equation(numbering: "(1)")
#set par(justify: true)
#set table(stroke: tablecolor)
#set heading(numbering: "1.1.")


#align(center, [#text(30pt)[
  Lo Básico] \ #text(13pt)[Apuntes Cap. 2 de Lyapunov Exponents - Arcady Pikovsky]
])
= Contexto matemático

Nos vamos a situar en la teoría de sistemas dinámicos en $RR^n$, cosa de que no sean objetos muy raros pero suficiente para modelar casi todo lo estudiado en física. Además, consideraremos que un sistema dinámico decente en palabras simples tiene una evolución que se determina por el lugar y tiempo en el instante "actual". Osea en verdad ignoramos que dependa de la historia del sistema, que es equivalente en verdad a nuestra primera restricción, ya que si el sistema depende de la historia, entonces habitualmente habita en un espacio de dimension infinita. Una excepción pueden ser sistemas dinámicos a tiempo discreto, donde la memoria finita no implica dimensión infinita, pero mejor cortemos esta distracción. Trabajaremos con 2 tipos principales de sistemas dinámicos.
#definición[Sistema dinámico][Un sistema dinámico viene dado por los siguientes ingredientes:
+ Un conjunto de estados (o espacio de fases) $M subset RR^n$. Hasta donde sé, no perdemos mucho si consideramos siempre $M = RR^n$.
+ Un conjunto temporal $cal(T)$.
  - Si $cal(T) = NN$ o $ZZ$, es a tiempo discreto.
  - Si $cal(T) = RR^+$ o $RR$, es a tiempo continuo.
+ Una función $U: cal(T) arrow.r M$, junto a su condición inicial $U_0 = U(0)$, representando la trayectoria del sistema.
+ Una función de evolución $F: M^n times cal(T) arrow M^n$ (en caso de ser discreto) o $F: M^n times cal(T) arrow RR^n$ (en caso de ser continua) tal que:
  - Caso continuo: $U(t + 1) = F(U(t), t)$.
  - Caso discreto: $dot(U)(t) = F(U(t), t)$.
]

La gracia de lo que vamos a hacer es estudiar la estabilidad de estas trayectorias ante perturbaciones. Osea queremos ver que si cambiamos la condición inicial un poco ¿qué irá a pasar con la norma de la diferencia $norm(u(t)) := norm(U(t) - tilde(U)(t))$?. Lo que usaremos es que lo más importante para esta pregunta (si $F$ es razonable) es la parte lineal de cómo evoluciona la diferencia. Esto porque incluso si uno se acerca debido a no linealidades, habrá un punto donde la parte lineal igual dominaría. Por esto, vamos a considerar otro sistema dinámico que sólo estudia esta parte lineal (doy estas justificaciones porque a futuro muchas veces no necesariamente consideraremos este sistema dinámico auxiliar como algo "chico" ya que a veces conviene que esté normalizado, o que lo dejemos crecer, debido a facilidades algorítmicas). El punto es que estudiar la parte lineal es lo importante. Y eso justo es lo que hacemos al considerar el sistema dinámico dado por:

$
  u(t + 1) = (partial F)/(partial U)(U(t), t) u(t) =: J_(U_0)(t)u
$

Esto me permite considerar una matriz que es dada por la multiplicación de las matrices desde 0 a t y obtener:
$
  u(t) = H(U_0, t) u(0)
$
O bien en el caso continuo:

$
  dot(u)(t) = (partial F)/(partial U)(U(t), t) u(t) =: K_(U_0)(t)u
$
La exponencial de una transformación lineal dim finita sigue siendo lineal así que la parte lineal a su vez cumple que:
$
  u(t) = H(U_0, t) u(0)
$

Osea que en verdad se parece montón a la formulación discreta. Aquí, $u(t) = U(t) - tilde(U)(t)$ ($tilde(U)(t)$ es la trayectoria ligeramente perturbada). Muchas veces en verdad uno estudia el caso discreto pero tiene propiedades casi iguales en el continuo. Hay un par de excepciones importantes sí.

= Exponentes de Lyapunov para $n = 1$, tiempo discreto

Digamos que estamos simplemente en $RR$, y en la versión discreta. Aquí podemos utilizar elementos de teoría ergódica para definir un exponente que nos caracteriza cómo se porta $u$. La idea es hacer la media ergódica de $ln(norm(J_(U_0)(t)))$ (logaritmo porque podemos separar las multiplicaciones en sumas, luego el $+ ln(u_0)$ lo podemos sacar de lo que estudiamos). Si $F$ nos define un sistema dinámico ergódico (en el sentido de medida), entonces debido al teorema ergódico puntual (el de Birkhoff para los alumbrados) esa media ergódica va a converger a la media en el sentido de medida, $angle.l ln(norm(J_(U_0)(t))) angle.r $. Este número, que nos caracteriza casi seguramente el comportamiento de la norma de la diferencia, se conoce como el *Exponente de Lyapunov* y lo denotamos $lambda$.

