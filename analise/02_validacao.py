#!/usr/bin/env python3
"""
Etapa A — validar o simulador contra a expectativa analítica.

Um simulador que não reproduz o caso com fórmula conhecida não merece
confiança no caso sem fórmula. Esta é a única etapa do projeto em que sabemos
a resposta certa de antemão, e por isso é a única chance de conferir.

Uso:  .venv/bin/python analise/02_validacao.py
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ))
from modelo.moran import (prob_fixacao_analitica, prob_fixacao_simulada,
                          trajetorias)

FIGS = RAIZ / "analise/figuras"
FIGS.mkdir(parents=True, exist_ok=True)

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def limpar(ax):
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)

def secao(t):
    print(f"\n{'─'*74}\n{t}\n{'─'*74}")


# ══════════════════════════════════════════════════════ teste 1: caso neutro
secao("TESTE 1 — CASO NEUTRO (s = 0)")
print("Sem vantagem nenhuma, a chance de um clone tomar a população deve ser")
print("exatamente a fatia que ele já ocupa: i₀/N. É o caso mais simples e o")
print("que pega erro de implementação mais grosseiro.\n")

print(f"{'N':>6}  {'i₀':>4}  {'esperado':>10}  {'simulado':>10}  {'±':>8}  {'desvios':>8}")
falhas = 0
for N, i0 in [(100, 1), (100, 10), (100, 50), (500, 1), (500, 25)]:
    esp = prob_fixacao_analitica(N, 0.0, i0)
    sim, err, _ = prob_fixacao_simulada(N, 0.0, i0, replicas=10_000, semente=42)
    z = abs(sim - esp) / err if err > 0 else 0
    marca = "ok" if z < 3 else "FALHOU"
    falhas += z >= 3
    print(f"{N:>6}  {i0:>4}  {esp:>10.5f}  {sim:>10.5f}  {err:>8.5f}  {z:>6.2f}σ  {marca}")

print("\n'desvios' é a distância entre simulado e esperado medida em erros")
print("padrão. Abaixo de 3σ é indistinguível de ruído amostral — o que se")
print("espera quando a implementação está certa.")


# ═══════════════════════════════════════════════════ teste 2: com vantagem
secao("TESTE 2 — COM VANTAGEM SELETIVA")
print("Agora a fórmula completa:  P_fix = (1 - (1/r)^i₀) / (1 - (1/r)^N),")
print("com r = 1+s. É aqui que seleção e deriva disputam de verdade.\n")

N = 200
print(f"N = {N}, partindo de uma única célula mutada\n")
print(f"{'s':>7}  {'esperado':>10}  {'simulado':>10}  {'±':>8}  {'desvios':>8}")
esses = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
sim_vals, err_vals = [], []
for s in esses:
    esp = prob_fixacao_analitica(N, s, 1)
    sim, err, _ = prob_fixacao_simulada(N, s, 1, replicas=10_000, semente=7)
    sim_vals.append(sim); err_vals.append(err)
    z = abs(sim - esp) / err if err > 0 else 0
    marca = "ok" if z < 3 else "FALHOU"
    falhas += z >= 3
    print(f"{s:>7.2f}  {esp:>10.5f}  {sim:>10.5f}  {err:>8.5f}  {z:>6.2f}σ  {marca}")


# ══════════════════════════════════════ teste 3: a regra de bolso, conferida
secao("TESTE 3 — A APROXIMAÇÃO PARA POPULAÇÃO GRANDE")
print("Para N grande e s pequeno, a fórmula colapsa em algo memorável.")
print("Convém conferir qual é o valor certo, porque circulam dois.\n")

N = 100_000
print(f"N = {N:,}\n")
print(f"{'s':>7}  {'exato':>10}  {'s/(1+s)':>10}  {'2s':>10}")
for s in [0.001, 0.005, 0.01, 0.05, 0.1]:
    exato = prob_fixacao_analitica(N, s, 1)
    print(f"{s:>7.3f}  {exato:>10.6f}  {s/(1+s):>10.6f}  {2*s:>10.6f}")

print("\nA aproximação correta para o processo de MORAN é s/(1+s) ≈ s.")
print("A regra '2s' que aparece em muito texto vem do modelo de")
print("Wright-Fisher, que tem variância de prole diferente. Os dois modelos")
print("descrevem a mesma biologia e discordam por um fator 2 — usar a")
print("aproximação errada dobraria nossa estimativa de aptidão.")


# ═════════════════════════════════════════════════════════════════ figuras
secao("FIGURAS")

# Fig 1 — simulado contra analítico. Duas séries: legenda presente.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ss = np.linspace(0, 0.55, 200)
ax.plot(ss, [prob_fixacao_analitica(200, s, 1) for s in ss],
        color=S1, linewidth=2, label="fórmula exata")
ax.errorbar(esses, sim_vals, yerr=np.array(err_vals)*1.96, fmt="o",
            color=S2, markersize=7, capsize=4, linewidth=1.6,
            label="simulação (10 mil réplicas, IC 95%)")
ax.set_xlabel("coeficiente de seleção  s")
ax.set_ylabel("probabilidade de fixação")
ax.set_title("O simulador reproduz a fórmula exata")
ax.text(0.97, 0.06, "N = 200, partindo de 1 célula mutada",
        transform=ax.transAxes, ha="right", fontsize=9, color=INK2)
ax.legend(loc="upper left")
limpar(ax)
fig.tight_layout(); fig.savefig(FIGS/"05_validacao_fixacao.png"); plt.close(fig)
print("  05_validacao_fixacao.png")

# Fig 2 — trajetórias: por que uma simulação não diz nada.
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
for ax, s, titulo in [(axes[0], 0.0, "sem vantagem (s = 0)"),
                      (axes[1], 0.1, "com vantagem (s = 0,1)")]:
    N = 200
    t, m = trajetorias(N, s, i0=1, n=40, passos=150_000, semente=1000)
    fixou = int((m[:, -1] >= N).sum())
    for k in range(m.shape[0]):
        venceu = m[k, -1] >= N
        ax.plot(t, m[k]/N, color=(S2 if venceu else S1),
                alpha=(0.85 if venceu else 0.30),
                linewidth=(1.4 if venceu else 1.0))
    ax.set_title(f"{titulo}   —   {fixou} de 40 fixaram", fontsize=10)
    ax.set_xlabel("tempo (gerações)")
    limpar(ax)
axes[0].set_ylabel("fração da população")
fig.suptitle("Mesmos parâmetros, destinos diferentes: a deriva decide",
             x=0.012, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIGS/"06_trajetorias.png"); plt.close(fig)
print("  06_trajetorias.png")


# ═════════════════════════════════════════════════════════════════ veredito
secao("VEREDITO")
if falhas == 0:
    print("Todos os testes passaram dentro de 3σ.")
    print("\nO simulador reproduz a fórmula exata no caso neutro e no caso com")
    print("seleção, para vários tamanhos de população e vários pontos de")
    print("partida. Está apto a ser usado onde não há fórmula.")
else:
    print(f"ATENÇÃO: {falhas} teste(s) fora de 3σ. Não prosseguir sem investigar.")
    sys.exit(1)
