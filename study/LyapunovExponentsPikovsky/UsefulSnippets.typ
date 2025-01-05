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


#align(center, [#text(30pt)[Salpicados Útiles] \ #text(13pt)[Apuntes de Lyapunov Exponents - Arcady Pikovsky]
])

#bloque[En este archivo busco recopilar secciones útiles esparcidas entre el fin del capítulo 2 y el inicio del capítulo 5 del libro.]

= Métodos para sistemas a tiempo continuo

Aparte de los métodos genéricos para encontrar exponentes de lyapunov, existen además métodos específicos para sistemas a tiempo continuo, basados en una descomposición adecuada de las matrices asociadas a nuestro sistema.  No me queda claro si esto se traduce en una ventaja computacional, así que hasta aquí llega mi comentario.

= Errores cerca de degeneracienes

En el libro se estudia el caso de un anillo de osciladores de Rössler, que para ciertos parámetros presenta un único exponente de lyapunov positivo. Además, gracias a la simetría traslacional del sistema, la mayoría de los exponentes vienen de a pares *pendiente: demostrar esto*. Para este caso, el atractor caótico es simplemente una trayectoria sincrónica, luego se dice que los modos de Fourier diagonalizan la ecuación de evolución de las perturbaciones. Yo veo esto como el paso justo antes de la MSF que ya estudiamos. Como resultado, la dinámica de cada modo de fourier está dada por una edo de sólo 3 variables. Cada modo de fourier entonces tiene asociado 3 exponentes de lyapunov, ya que el iscilador de rossler tiene 3 variables. Súper. 
 

