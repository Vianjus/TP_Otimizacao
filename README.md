## Modelo AMPL para Gusek do [Artigo](./artigo_seminario.pdf)
---
### 1. Função Objetivo

**Artigo:**  
Minimiza o custo de instalação/manutenção dos bins + custo operacional dos veículos (tempo total de serviço dos veículos).

$$
\min \sum_{b \in B} (CIN_b \sum_{i \in I} n_{bi}) + CCV \sum_{t \in T} \sum_{v \in V} TT_{vt}
$$

**Meu modelo:**  
```mod
minimize Total_Cost:
    sum{b in B, i in I} CINb[b] * nb[b,i]
    + CCV * sum{t in T, v in V} (
        TU * sum{i in I} x[i,0,v,t]
        + sum{i in I union {0}, j in I union {0}} x[i,j,v,t] * C[i,j]
        + sum{i in I union {0}, j in I union {0}, b in B} Sb[b] * z[i,j,b,v,t]
      );
```
---

### 2. Variáveis de Decisão

**Artigo:**  
- $x_{ijvt}$: se veículo $v$ vai de $i$ para $j$ no dia $t$
- $y_{ijvt}$: carga transportada de $i$ para $j$ por $v$ no dia $t$
- $w_{it}$: resíduo acumulado no ponto $i$ no dia $t$
- $w_{i}^{max}$: máximo resíduo acumulado no ponto $i$
- $n_{bi}$: se combinação $b$ de bins está no ponto $i$
- $z_{ijbvt}$: linearização do tempo de serviço

**Meu modelo:**
```mod
# Routing variables
var x{i in I union {0}, j in I union {0}, v in V, t in T}, binary;
var y{i in I union {0}, j in I union {0}, v in V, t in T}, >= 0;

# Waste accumulation variables
var w{i in I, t in T}, >= 0;
var wmax{i in I}, >= 0;

# Bin allocation variables
var nb{b in B, i in I union {0}}, binary;

# Linearization variables for service time (Glover, 1975)
var z{i in I union {0}, j in I union {0}, b in B, v in V, t in T} >= 0;
```
---

### 3. Restrições

**Artigo:**  
- (2a) Não pode ter bin no depósito
- (2b) Exatamente uma combinação de bin por ponto
- (2c) Capacidade dos bins ≥ máximo resíduo acumulado
- (2d) Não pode ter self-loop
- (2e) Não coleta em dias de descanso
- (2f) Conservação de fluxo
- (2g) Máximo de uma rota por veículo por dia
- (2h) Tempo de serviço ≤ jornada
- (2i) Carga ≤ capacidade do veículo
- (2j) Balanceamento de carga/subtour elimination
- (2m) Resíduo ≤ capacidade máxima
- (3a)-(3c) Linearização da acumulação de resíduos (cíclica)
- (4a)-(4c) Linearização do tempo de serviço (Glover)
- (6a)-(6b) Quebra de simetria

**Meu modelo:**  
Todas essas restrições estão implementadas, com nomes e fórmulas equivalentes.

---

### 4. Linearização e Simetria

**Artigo:**  
- Linearização do tempo de serviço com variáveis $z_{ijbvt}$ (Eq. 4a-4c)
- Quebra de simetria (Eq. 6a-6b)

**Meu modelo:**  
Inclui as mesmas linearizações e cortes de simetria.

---

### 5. Parâmetros

**Artigo:**  
Todos os parâmetros do artigo estão presentes no meu modelo (capacidades, custos, tempos, etc.).

---
