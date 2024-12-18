#import "@preview/ctheorems:1.1.3": *
#import "@preview/physica:0.9.3": *
#import "@preview/cetz:0.3.1"
#import "@preview/cetz-plot:0.1.0"
#import "@preview/xarrow:0.3.1": xarrow

#let (backgroundcolor, headingcolors, textcolors, tablecolor) = (
  rgb("#282A36"),
  (rgb("#BD93F9"), rgb("#FF79C6"), rgb("#8BE9FD"), rgb("#50FA7B")),
  (rgb("#F8F8F2"), rgb("#F1FA8C"), rgb("#FFB86C")),
  rgb("#6272A4")
)

#let teorema = thmbox("teorema", "Teorema")
#let lema = thmbox("lema", "Lema")
#let corolario = thmplain(
  "corolario",
  "Corolario",
  base: "teorema"
)

#let definición = thmbox("definición", "Definición")
#let ejemplo = thmplain("ejemplo", "Ejemplo").with(numbering: none)
#let dem = thmproof("demostración", "Demostración")
#let propiedad = thmplain(
  "propiedad",
  "Propiedad",
  base: "definición"
)
#let obs = thmplain("obs", "Obs").with(numbering: none)

#let innerproduct(x, y) = {
  [$lr(angle.l #x, #y angle.r)$]
}

#let braket(a, b) = {
  [$lr(angle.l #a|#b angle.r)$]
}

#let bra(x) = {
  [$lr(angle.l #x|)$]
}

#let ket(x) = {
  [$lr(|#x angle.r)$]
}

#let implies = [$=>$]

#let iff = [$<=>$]

#let colMath(x, color) = text(fill: color)[$#x$]

#let bloque(BlockText, Titulo: none) = {[
  *#Titulo*
  #block(stroke: textcolors.at(2), inset: 1em)[#BlockText]
]
}