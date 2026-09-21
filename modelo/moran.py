"""
Processo de Moran com seleção — o modelo mínimo de competição clonal.

A REGRA, em cada passo de tempo:
  1. uma célula é escolhida para se reproduzir, com probabilidade
     proporcional à aptidão (mutantes pesam r = 1+s, normais pesam 1)
  2. uma célula é escolhida para morrer, uniformemente ao acaso
  3. a população volta ao tamanho original

O passo 1 é a seleção. O passo 2 é a deriva. As duas forças emergem da mesma
regra, sem serem programadas separadamente.

DUAS DECISÕES DE IMPLEMENTAÇÃO QUE MUDAM TUDO DE VELOCIDADE

1. Não rastreamos células. Como só existem dois tipos, o estado do sistema é
   um único número: `i`, quantas células são mutantes. O modelo vira uma
   cadeia de nascimento-e-morte — exata, e sem custo proporcional a N.

2. Não rodamos uma réplica por vez. Todas as réplicas avançam juntas num
   array do numpy. O número de operações em Python passa a ser o número de
   PASSOS, não passos × réplicas. Na prática isso vale duas ordens de
   grandeza.
"""

import numpy as np


def prob_fixacao_analitica(N, s, i0=1):
    """
    Probabilidade de o clone tomar a população inteira, em forma fechada.

    Resultado clássico do processo de Moran:

                    1 - (1/r)^i0
        P_fix  =  ────────────────       com r = 1 + s
                    1 - (1/r)^N

    No caso neutro (s = 0) o limite é simplesmente i0/N: sem vantagem, a
    chance de vencer é a fatia que você já ocupa.

    É esta fórmula que o simulador precisa reproduzir.
    """
    if s == 0:
        return i0 / N
    x = 1.0 / (1.0 + s)
    return (1.0 - x**i0) / (1.0 - x**N)


def _passo_vetorizado(i, N, s, rng):
    """
    Avança um passo da cadeia mergulhada, para todas as réplicas ativas.

    A CADEIA MERGULHADA
    A maioria dos passos de Moran não muda nada — um mutante substitui outro
    mutante. Em vez de sortear esses passos inúteis, sorteamos só a DIREÇÃO
    da próxima mudança, condicionada a que uma mudança ocorra:

        P(subir | mudou) = p_sobe / (p_sobe + p_desce)

    Distribuição idêntica, muito mais rápido. O preço é perder a noção de
    tempo real — irrelevante quando a pergunta é "fixa ou não".
    """
    r = 1.0 + s
    x = i.astype(np.float64)
    peso = r * x + (N - x)
    p_sobe = (r * x / peso) * ((N - x) / N)
    p_desce = ((N - x) / peso) * (x / N)
    p = p_sobe / (p_sobe + p_desce)
    return np.where(rng.random(x.size) < p, 1, -1)


def prob_fixacao_simulada(N, s, i0=1, replicas=10_000, semente=None,
                          max_passos=50_000_000):
    """
    Estima a probabilidade de fixação rodando muitas realizações em paralelo.

    Devolve (estimativa, erro_padrao, passos_ate_absorver_todas).

    O erro padrão vem da binomial — sqrt(p(1-p)/n). É quanto a estimativa
    oscila só por causa do número finito de réplicas, e serve para decidir se
    uma divergência da fórmula é real ou é ruído amostral.
    """
    rng = np.random.default_rng(semente)
    i = np.full(replicas, i0, dtype=np.int64)
    ativo = (i > 0) & (i < N)
    passos = 0

    while ativo.any():
        idx = np.flatnonzero(ativo)
        i[idx] += _passo_vetorizado(i[idx], N, s, rng)
        ativo = (i > 0) & (i < N)
        passos += 1
        if passos > max_passos:
            raise RuntimeError(f"não absorveu em {passos} passos (N={N}, s={s})")

    p = (i >= N).mean()
    return p, np.sqrt(max(p * (1 - p), 1e-12) / replicas), passos


def trajetorias(N, s, i0=1, n=40, passos=200_000, semente=None):
    """
    Roda `n` realizações com o tempo explícito, guardando o caminho de cada.

    Diferente da versão acima, aqui cada passo é um passo de Moran de verdade,
    incluindo os que não mudam nada — porque para desenhar a trajetória o
    tempo importa.

    Devolve (tempos, matriz) com matriz de forma (n, len(tempos)). Réplicas
    já absorvidas ficam congeladas no valor final, que é o comportamento
    correto: extinto continua extinto.

    O tempo está em gerações, e N passos de Moran ≈ uma geração.
    """
    rng = np.random.default_rng(semente)
    r = 1.0 + s
    i = np.full(n, i0, dtype=np.int64)
    hist = [i.copy()]

    for _ in range(passos):
        ativo = (i > 0) & (i < N)
        if not ativo.any():
            break
        idx = np.flatnonzero(ativo)
        x = i[idx].astype(np.float64)
        peso = r * x + (N - x)
        p_sobe = (r * x / peso) * ((N - x) / N)
        p_desce = ((N - x) / peso) * (x / N)
        u = rng.random(x.size)
        delta = np.where(u < p_sobe, 1, np.where(u < p_sobe + p_desce, -1, 0))
        i[idx] += delta
        hist.append(i.copy())

    m = np.array(hist).T                      # (n, passos)
    t = np.arange(m.shape[1]) / N
    return t, m
