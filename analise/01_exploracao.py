#!/usr/bin/env python3
"""
Exploração do conjunto agregado de Watson et al. 2020.

Objetivo: entender o que o dado permite responder ANTES de escrever o modelo.
Um simulador calibrado contra uma observável que o dado não mede é trabalho
perdido.

Uso:  .venv/bin/python analise/01_exploracao.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RAIZ = Path(__file__).resolve().parent.parent
CSV = (RAIZ / "referencia/dados/watson2020/Maximum_likelihood_estimations"
            / "Maximum likelihood estimations - data files"
            / "all_studies_trimmed_all_genes.csv")
FIGS = RAIZ / "analise/figuras"
FIGS.mkdir(parents=True, exist_ok=True)

# paleta: superfície clara, tinta recessiva, três séries no máximo em dispersão
SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "axes.labelsize": 10,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def limpar_eixos(ax, grid="y"):
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(axis=grid, which="major")
    ax.grid(axis=("x" if grid == "y" else "y"), which="major", visible=False)

def secao(t):
    print(f"\n{'─'*72}\n{t}\n{'─'*72}")


# ──────────────────────────────────────────────────────────── carga e limpeza
secao("1. CARGA E LIMPEZA")

bruto = pd.read_csv(CSV, lineterminator="\r")
bruto.columns = [c.strip() for c in bruto.columns]
print(f"linhas lidas: {len(bruto)}")

# a idade vem como texto: 663 linhas trazem 'noagedata'. Converter com
# errors='coerce' transforma o inválido em NaN em vez de explodir.
bruto["age"] = pd.to_numeric(bruto["age"], errors="coerce")
bruto["VAF"] = pd.to_numeric(bruto["VAF"], errors="coerce")

df = bruto.dropna(subset=["age", "VAF"]).copy()
print(f"descartadas por falta de idade: {len(bruto) - len(df)}")
print(f"utilizáveis: {len(df)}")
print(f"\nidade:  {df.age.min():.0f} a {df.age.max():.0f} anos "
      f"(mediana {df.age.median():.0f})")
print(f"VAF:    {df.VAF.min():.4f} a {df.VAF.max():.4f} "
      f"(mediana {df.VAF.median():.4f})")

# quem ficou de fora importa: se as linhas sem idade forem de um estudo só,
# a perda é sistemática, não aleatória.
perdidas = bruto[bruto.age.isna()]
if len(perdidas):
    print("\ncoortes perdidas por falta de idade:")
    for est, n in perdidas.study.value_counts().items():
        tot = (bruto.study == est).sum()
        print(f"  {est:<16} {n:>4} de {tot:>4}  ({100*n/tot:.0f}%)")


# ─────────────────────────────────────────── o viés que precisa ser encarado
secao("2. LIMITE DE DETECÇÃO POR COORTE")
print("Cada estudo sequenciou com profundidade diferente, então cada um enxerga")
print("clones a partir de um tamanho mínimo diferente. Comparar VAF entre")
print("coortes sem levar isso em conta mistura biologia com instrumento.\n")

lim = (df.groupby("study")
         .agg(n=("VAF", "size"), vaf_min=("VAF", "min"),
              vaf_p05=("VAF", lambda s: s.quantile(0.05)),
              vaf_mediana=("VAF", "median"), idade_med=("age", "median"))
         .sort_values("vaf_min"))
print(lim.to_string(float_format=lambda x: f"{x:.4f}"))


# ──────────────────────────────────────────────────── a relação central
secao("3. VAF CONTRA IDADE")
print("Sob crescimento exponencial de clone, VAF ~ exp(s·t): log(VAF) deve")
print("crescer linearmente com a idade, e a inclinação carrega a aptidão s.\n")

x, y = df.age.values, np.log(df.VAF.values)
inc, interc, r, p, err = stats.linregress(x, y)
print(f"regressão de log(VAF) sobre idade, todos os genes:")
print(f"  inclinação: {inc:+.5f} por ano  (erro padrão {err:.5f})")
print(f"  p-valor:    {p:.3e}")
print(f"  R²:         {r**2:.4f}")
print(f"\nR² baixo é esperado e informativo: a idade sozinha não determina o")
print(f"tamanho do clone. Cada variante tem sua própria aptidão, e o acaso da")
print(f"deriva domina clones pequenos. É justamente o que um modelo estocástico")
print(f"precisa reproduzir — dispersão, não uma curva.")

print("\ninclinação por gene (apenas genes com n ≥ 40):")
linhas = []
for gene, g in df.groupby("gene"):
    if len(g) < 40:
        continue
    s, _, rr, pp, ee = stats.linregress(g.age.values, np.log(g.VAF.values))
    linhas.append((gene, len(g), s, ee, pp, rr**2))
for gene, n, s, ee, pp, r2 in sorted(linhas, key=lambda t: -t[2]):
    sig = "*" if pp < 0.05 else " "
    print(f"  {gene:<8} n={n:>4}  inclinação {s:+.5f} ± {ee:.5f}{sig}  R²={r2:.3f}")
print("  (* p < 0,05)")


# ──────────────────────────────────────────── a pergunta biológica concreta
secao("4. DNMT3A R882 CONTRA O RESTO DE DNMT3A")
print("R882 é o ponto quente mais conhecido de CHIP. Se ele confere vantagem")
print("maior, os clones devem ser maiores na mesma idade.\n")

d3 = df[df.gene == "DNMT3A"].copy()
d3["grupo"] = np.where(d3.variant.str.startswith("R882"), "R882", "outras")
for grp, g in d3.groupby("grupo"):
    print(f"  {grp:<8} n={len(g):>4}  VAF mediana {g.VAF.median():.4f}  "
          f"idade mediana {g.age.median():.0f}")

a = d3[d3.grupo == "R882"].VAF.values
b = d3[d3.grupo == "outras"].VAF.values
u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
print(f"\n  Mann-Whitney U: p = {pu:.3e}")
print(f"  (teste não paramétrico: as distribuições de VAF são assimétricas,")
print(f"   então comparar médias seria enganoso)")


# ────────────────────────────────────────────────────────────────── figuras
secao("5. FIGURAS")

# Fig 1 — a relação central. Série única: sem legenda, o título nomeia.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.scatter(df.age, df.VAF, s=14, alpha=0.45, color=S1,
           edgecolors="none", rasterized=True)
xx = np.linspace(df.age.min(), df.age.max(), 100)
ax.plot(xx, np.exp(interc + inc*xx), color=INK, linewidth=2,
        label="ajuste exponencial")
ax.set_yscale("log")
ax.set_xlabel("idade (anos)")
ax.set_ylabel("frequência alélica (VAF)")
ax.set_title("Clones detectados crescem com a idade, mas a dispersão domina")
ax.text(0.98, 0.04, f"n = {len(df)} variantes\ninclinação {inc:+.4f}/ano   R² = {r**2:.3f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=INK2)
limpar_eixos(ax, grid="both")
ax.grid(axis="x", visible=True)
fig.tight_layout(); fig.savefig(FIGS/"01_vaf_por_idade.png"); plt.close(fig)
print("  01_vaf_por_idade.png")

# Fig 2 — o viés instrumental. Série única, ordenada por magnitude.
fig, ax = plt.subplots(figsize=(7.2, 4.0))
o = lim.sort_values("vaf_min")
pos = np.arange(len(o))
ax.hlines(pos, o.vaf_min, o.vaf_mediana, color=BASE, linewidth=2)
ax.scatter(o.vaf_min, pos, s=64, color=S1, zorder=3, label="menor VAF detectada")
ax.scatter(o.vaf_mediana, pos, s=64, color=S2, zorder=3, label="VAF mediana")
ax.set_xscale("log")
ax.set_yticks(pos); ax.set_yticklabels([f"{i}  (n={int(n)})" for i, n in zip(o.index, o.n)])
ax.set_xlabel("frequência alélica (VAF, escala log)")
ax.set_title("Cada coorte enxerga a partir de um tamanho de clone diferente")
ax.legend(loc="lower right")
limpar_eixos(ax, grid="x")
fig.tight_layout(); fig.savefig(FIGS/"02_limite_deteccao.png"); plt.close(fig)
print("  02_limite_deteccao.png")

# Fig 3 — duas séries: legenda presente e rótulo direto.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for grp, cor in (("R882", S2), ("outras", S1)):
    v = np.sort(d3[d3.grupo == grp].VAF.values)
    ax.step(v, np.arange(1, len(v)+1)/len(v), where="post",
            color=cor, linewidth=2, label=f"DNMT3A {grp} (n={len(v)})")
ax.set_xscale("log")
ax.set_xlabel("frequência alélica (VAF, escala log)")
ax.set_ylabel("proporção acumulada")
ax.set_title("R882 não produz clones maiores neste conjunto" if pu > 0.05
             else "R882 produz clones maiores que as demais variantes de DNMT3A")
ax.text(0.03, 0.95, f"Mann-Whitney p = {pu:.1e}", transform=ax.transAxes,
        va="top", fontsize=9, color=INK2)
ax.legend(loc="lower right")
limpar_eixos(ax, grid="both")
fig.tight_layout(); fig.savefig(FIGS/"03_dnmt3a_r882.png"); plt.close(fig)
print("  03_dnmt3a_r882.png")

# Fig 4 — pequenos múltiplos: evita pintar 9 coortes com 9 cores.
top = df.gene.value_counts().head(6).index.tolist()
fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.0), sharex=True, sharey=True)
for ax, gene in zip(axes.ravel(), top):
    g = df[df.gene == gene]
    ax.scatter(g.age, g.VAF, s=12, alpha=0.5, color=S1, edgecolors="none")
    if len(g) >= 10:
        s, i2, rr, pp, _ = stats.linregress(g.age.values, np.log(g.VAF.values))
        xs = np.linspace(g.age.min(), g.age.max(), 50)
        ax.plot(xs, np.exp(i2 + s*xs), color=INK, linewidth=1.6)
        ax.set_title(f"{gene}   n={len(g)}   {s:+.4f}/ano", fontsize=10)
    ax.set_yscale("log")
    limpar_eixos(ax, grid="both")
for ax in axes[1]:
    ax.set_xlabel("idade (anos)")
for ax in axes[:, 0]:
    ax.set_ylabel("VAF")
fig.suptitle("Cada gene condutor tem sua própria taxa de crescimento",
             x=0.012, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGS/"04_por_gene.png"); plt.close(fig)
print("  04_por_gene.png")


# ────────────────────────────────────────────────────── o que o dado permite
secao("6. O QUE ESTE DADO PERMITE — E O QUE NÃO PERMITE")
print("""PERMITE
  · distribuição de VAF condicionada à detecção, por idade e por gene
  · comparar taxas de crescimento entre genes condutores
  · estimar aptidão a partir da inclinação de log(VAF) sobre idade

NÃO PERMITE
  · prevalência de CHIP na população: o arquivo traz apenas as variantes
    DETECTADAS, sem o denominador de quantas pessoas foram rastreadas por
    faixa de idade. Prevalência exige esse denominador, e ele não está aqui.
  · trajetória de um mesmo indivíduo ao longo do tempo: cada linha é uma
    medição única. A dinâmica temporal é inferida comparando pessoas de
    idades diferentes, não acompanhando as mesmas.

CONSEQUÊNCIA PARA O MODELO
  A observável a reproduzir é a DISTRIBUIÇÃO de VAF por idade, truncada pelo
  limite de detecção de cada coorte — não uma curva média. O simulador tem de
  gerar dispersão realista, e a comparação com o observado precisa aplicar o
  mesmo truncamento, senão compara maçã com laranja.""")
