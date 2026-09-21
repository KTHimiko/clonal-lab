# Referências — o que existe e como obter

Levantamento do que está publicamente disponível sobre dinâmica clonal em
hematopoiese, feito antes de escrever qualquer modelo. Conhecer o formato do
dado real evita construir um simulador que não conversa com nada.

> [!IMPORTANT]
> Nada aqui é código nosso. Este diretório é material de terceiros, para
> leitura e referência, e **não deve ser versionado** junto com o projeto.

---

## O que foi obtido

### Mitchell et al. 2022 — código

**Clonal dynamics of haematopoiesis across the human lifespan**, *Nature*
607:343–350. Mapeia a dinâmica clonal ao longo da vida humana e mostra a
mudança de regime após os ~70 anos, com poucos clones dominando.

- Repositório: <https://github.com/emily-mitchell/normal_haematopoiesis>
- Clonado em `codigo/normal_haematopoiesis` — 27 MB, R
- Artigo aberto: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9177428/>

> [!CAUTION]
> **O repositório não declara licença.** Sem licença, o padrão legal é "todos
> os direitos reservados": pode-se ler e aprender, **não** reutilizar o código
> num projeto próprio. Se algum trecho for útil, reimplemente a partir do
> método descrito no artigo, não por cópia.

**O que ele ensina, e é o mais valioso:**

| Diretório | Conteúdo |
|---|---|
| `6_population_modelling` | estimativa do tamanho populacional de células-tronco |
| `7_phylofit` | ajuste de dinâmica a partir de filogenias |
| `8_driver_modelling` | **ABC para inferir efeito de aptidão de mutações condutoras** |
| `10_simulations_for_figures` | simulações que geram as figuras |
| `11_LOY_simulations` | perda do cromossomo Y |

A metodologia central é **ABC** (Computação Bayesiana Aproximada): simular
milhares de vezes variando parâmetros, calcular estatísticas-resumo de cada
execução, e manter as que chegam perto do observado. A distribuição dos
parâmetros mantidos é a estimativa.

É massivamente paralela — cada simulação independe das outras — e portanto o
tipo exato de carga para a qual um cluster existe.

Eles submetem com `bsub`, o LSF do Sanger. Nós usamos SLURM. A metodologia é
agnóstica ao escalonador; muda só a linha de submissão.

O motor de simulação deles é o pacote `Rsimpop`.

---

## O que exige download manual

Estas fontes exigem navegador. Tentativas programáticas esbarram em
autenticação ou proteção contra robô — documentado aqui para não repetir o
esforço.

### Watson et al. 2020 — dados

**The evolutionary dynamics and fitness landscape of clonal hematopoiesis**,
*Science* 367:1449–1454. Estimou o efeito de aptidão de variantes condutoras a
partir de dados de ~50.000 pessoas, e delimitou o número de células-tronco
hematopoiéticas.

- Dados: <https://doi.org/10.5061/dryad.83bk3j9mw> — **licença CC0**
  (domínio público, uso irrestrito)
- Arquivo: `Age_prevalence_of_DNMT3A_R882H_and_R882C_variants.zip`, 57 KB
- Artigo (PDF aberto):
  <https://web.stanford.edu/group/dsfisher/papers/pdf/watson_et_al_2020.pdf>

*A API do Dryad responde `Unauthorized, must have current bearer token` e o
endpoint direto devolve 403. Baixe pela página.*

Salve em `dados/watson2020/`.

### Mitchell et al. 2022 — matrizes de dados

- Mendeley Data: <https://data.mendeley.com/datasets/np54zjkvxr/1>
- Dados brutos de sequenciamento: EGA, acesso **controlado** (EGAD00001007851)
  — fora do nosso escopo, exige aprovação institucional

Salve em `dados/mitchell2022/`.

---

## Leitura de base

Ordem sugerida, do fenômeno ao método:

1. **Jaiswal et al. 2014**, *NEJM* — estabelece CHIP como fenômeno e sua
   associação com mortalidade e risco cardiovascular. É a porta de doença do
   nicho.
2. **Genovese et al. 2014**, *NEJM* — mesma descoberta, independente.
3. **Watson et al. 2020**, *Science* — quantifica aptidão a partir de
   frequência alélica. É o alvo de calibração do nosso modelo.
4. **Mitchell et al. 2022**, *Nature* — dinâmica ao longo da vida, e a
   metodologia ABC que queremos reproduzir.

---

## O que isso define para o projeto

**A observável de calibração** é a distribuição de frequência alélica
(*variant allele frequency*) em função da idade. É o que os dados públicos
trazem e o que o modelo precisa prever.

**O método de inferência** é ABC: não há verossimilhança fechada para um
modelo baseado em agentes, então compara-se por estatísticas-resumo.

**A carga computacional** é a varredura: milhares de simulações independentes.
É onde o cluster entra, e é a razão de ele existir neste projeto.

**A restrição de dados** se confirma: tudo que precisamos é aberto ou CC0. O
único material controlado (EGA) é sequenciamento bruto, que não usaríamos de
qualquer forma.
