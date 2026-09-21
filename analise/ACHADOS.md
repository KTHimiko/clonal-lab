# Achados da exploração — Watson et al. 2020

Reproduzível com `.venv/bin/python analise/01_exploracao.py`.
Figuras em `analise/figuras/`.

---

## 1. A perda de dados é estrutural, não aleatória

663 das 1.674 linhas não têm idade. Mas não é uma perda difusa:

| Coorte | Perdidas | Total | |
|---|---|---|---|
| Jaiswal2014 | 466 | 466 | **100%** |
| Genovese2014 | 197 | 197 | **100%** |

**As duas maiores coortes foram excluídas por inteiro.** Nenhuma outra perdeu
uma linha sequer.

Isso muda a interpretação: não estamos com uma amostra 40% menor do mesmo
universo — estamos com um universo diferente, sem os dois estudos
populacionais clássicos de CHIP. Restam **1.012 variantes**, e qualquer
conclusão vale para as sete coortes que sobraram, não para o conjunto que o
artigo analisou.

> Se algum dia precisarmos de Jaiswal e Genovese, a idade está nos artigos
> originais — só não nesta tabela agregada.

---

## 2. O limite de detecção varia 300 vezes entre coortes

| Coorte | n | menor VAF | VAF mediana | idade mediana |
|---|---|---|---|---|
| Young2019 | 158 | 0,0008 | 0,0019 | 61 |
| Young2016 | 26 | 0,0011 | 0,0034 | 58 |
| Acuna2017 | 182 | 0,0016 | 0,0066 | 57 |
| McKerrel2015 | 112 | 0,0079 | 0,0272 | 78 |
| Desai2018 | 31 | 0,0182 | 0,0344 | 70 |
| Coombs2017 | 378 | 0,0192 | 0,0509 | 70 |
| **ZinkWGS** | 125 | **0,2400** | 0,3300 | 71 |

O ZinkWGS é sequenciamento de genoma inteiro — profundidade baixa, então só
enxerga clones que já tomaram um quarto do sangue. O Young2019 usa
sequenciamento dirigido e profundo, e vê clones 300 vezes menores.

**Consequência direta para o modelo:** comparar a distribuição simulada com a
observada exige aplicar, na simulação, o mesmo truncamento de cada coorte.
Sem isso, o modelo é penalizado por prever clones pequenos que o instrumento
nunca poderia ter visto.

---

## 3. O sinal de crescimento é real, e ruidoso do jeito certo

Regressão de log(VAF) sobre idade, todas as 1.012 variantes:

```
inclinação  +0,0387 por ano   (erro padrão 0,0039)
p-valor      3,4 × 10⁻²²
R²           0,089
```

Altamente significativo e com poder explicativo baixíssimo. Isso **não** é
contradição — é a assinatura de um processo estocástico.

A idade não determina o tamanho do clone: determina a distribuição de
tamanhos possíveis. Cada variante tem sua própria aptidão, e a deriva domina
enquanto o clone é pequeno. Um modelo determinístico ajustaria a linha e
erraria tudo o que importa; **o que precisamos reproduzir é a nuvem, não a
reta.**

Esse é o argumento mais forte a favor de simulação baseada em agentes em vez
de fórmula fechada.

---

## 4. Cada gene condutor cresce a uma taxa própria

| Gene | n | inclinação (por ano) | R² | |
|---|---|---|---|---|
| TET2 | 95 | **+0,0969 ± 0,0141** | 0,338 | significativo |
| DNMT3A | 421 | **+0,0468 ± 0,0063** | 0,116 | significativo |
| JAK2 | 45 | +0,0195 ± 0,0110 | 0,068 | não significativo |

TET2 cresce cerca de **duas vezes mais rápido** que DNMT3A nesta amostra, e
com R² três vezes maior — a relação idade-tamanho é bem mais limpa nele.

JAK2 não atinge significância, mas com n=45 isso diz mais sobre o tamanho da
amostra que sobre a biologia.

Isto define um requisito do modelo: **aptidão é por variante, não global.**
Um simulador com um único `s` não consegue reproduzir três genes com taxas
diferentes.

---

## 5. R882 confere vantagem mensurável

O ponto quente mais conhecido de CHIP, contra as demais variantes do mesmo
gene:

| Grupo | n | VAF mediana | idade mediana |
|---|---|---|---|
| DNMT3A R882 | 96 | **0,0390** | 63 |
| DNMT3A outras | 325 | 0,0215 | 63 |

Mann-Whitney U: **p = 5,4 × 10⁻⁵**

**Mesma idade mediana, clones com quase o dobro do tamanho.** Como a idade
está controlada por construção, a diferença é atribuível à variante, não ao
tempo de exposição.

É um alvo de validação excelente: um modelo bem calibrado deve inferir um `s`
maior para R882 que para o resto de DNMT3A, e a razão entre os dois deve ser
comparável ao que se observa aqui.

---

## 6. O que este dado não permite

**Prevalência de CHIP na população.** O arquivo traz apenas as variantes
detectadas, sem o denominador de quantas pessoas foram rastreadas por faixa
etária. Prevalência exige esse denominador, e ele não está aqui — nem pode ser
inferido, porque cada coorte rastreou populações de tamanhos diferentes.

**Trajetória individual.** Cada linha é uma medição única de uma pessoa. A
dinâmica temporal é inferida comparando pessoas de idades diferentes, não
acompanhando as mesmas ao longo do tempo. Isso é um estudo transversal, e
carrega todas as limitações de um: efeito de coorte, viés de sobrevivência, e
a impossibilidade de observar clones que desapareceram.

---

## O que isso decide sobre o modelo

| Achado | Requisito que impõe |
|---|---|
| Dispersão domina (R² = 0,09) | modelo estocástico, não determinístico |
| Limite de detecção varia 300× | truncar a simulação por coorte antes de comparar |
| Genes crescem a taxas diferentes | aptidão por variante, não um `s` global |
| R882 tem quase o dobro do tamanho | alvo de validação quantitativo |
| Sem denominador populacional | calibrar por distribuição de VAF, nunca por prevalência |
| Estudo transversal | comparar populações simuladas por idade, não trajetórias |

**A observável de calibração, definida com precisão:** a distribuição de VAF
entre clones detectados, estratificada por idade e por gene, truncada no
limite de detecção da coorte correspondente.
