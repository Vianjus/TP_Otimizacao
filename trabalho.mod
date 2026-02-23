# ============================================================================
# Modelagem e resolucao integrada de roteamento periodico de veiculos
# e localizacao de instalacoes com capacidade no contexto de
# coleta de residuos solidos
# ============================================================================

set I;                          # Pontos de coleta (excluindo o deposito)
set B;                          # Combinacoes de contentores
set V;                          # Veiculos
set T;                          # Dias no horizonte de planejamento
set Tprime within T;            # Dias de descanso (sem coleta)

param Q;                        # Capacidade do veiculo (m3)
param C{i in I union {0}, j in I union {0}};  # Tempo de viagem entre pontos (min)
param Sb{b in B};                # Tempo de servico da combinacao b (min)
param TU;                        # Tempo de descarga no deposito (min)
param Wi{i in I};                 # Geracao diaria de residuos no ponto i (m3/dia)
param CAPb{b in B};               # Capacidade da combinacao b de contentores (m3)
param CINb{b in B};               # Custo de instalacao e manutencao da combinacao b (US$)
param CCV;                        # Custo operacional dos veiculos (US$/min)
param TL;                         # Duracao da jornada de trabalho (min)
param BigM := max{b in B} CAPb[b]; # Big M para linearizacao do acumulo de residuos

# ============================================================================
# VARIAVEIS DE DECISAO
# ============================================================================

# Variaveis de roteamento
var x{i in I union {0}, j in I union {0}, v in V, t in T}, binary;
var y{i in I union {0}, j in I union {0}, v in V, t in T}, >= 0;

# Variaveis de acumulo de residuos
var w{i in I, t in T}, >= 0;
var wmax{i in I}, >= 0;

# Variaveis de alocacao de contentores
var nb{b in B, i in I union {0}}, binary;

# Variaveis de linearizacao para tempo de servico (Glover, 1975)
var z{i in I union {0}, j in I union {0}, b in B, v in V, t in T} >= 0;

# ============================================================================
# FUNCAO OBJETIVO
# ============================================================================

minimize Total_Cost:
    sum{b in B, i in I} CINb[b] * nb[b,i]
    + CCV * sum{t in T, v in V} (
        TU * sum{i in I} x[i,0,v,t]
        + sum{i in I union {0}, j in I union {0}} x[i,j,v,t] * C[i,j]
        + sum{i in I union {0}, j in I union {0}, b in B} Sb[b] * z[i,j,b,v,t]
      );

# ============================================================================
# RESTRICOES
# ============================================================================

# (2a) Nenhuma combinacao de contentor no deposito
s.t. c2a{b in B}:
    nb[b,0] = 0;

# (2b) Exatamente uma combinacao de contentor por ponto de coleta
s.t. c2b{i in I}:
    sum{b in B} nb[b,i] = 1;

# (2c) Capacidade dos contentores >= maximo residuo acumulado
s.t. c2c{i in I}:
    sum{b in B} CAPb[b] * nb[b,i] >= wmax[i];

# (2d) Sem auto-loops
s.t. c2d{i in I union {0}, v in V, t in T}:
    x[i,i,v,t] = 0;

# (2e) Sem coleta nos dias de descanso
s.t. c2e{v in V, t in Tprime, i in I union {0}, j in I union {0}}:
    x[i,j,v,t] = 0;

# (2f) Conservacao de fluxo
s.t. c2f{j in I union {0}, v in V, t in T}:
    sum{i in I union {0}} x[i,j,v,t] = sum{i in I union {0}} x[j,i,v,t];

# (2g) No maximo uma rota por veiculo por dia (exceto dias de descanso)
s.t. c2g{v in V, t in T diff Tprime}:
    sum{i in I} x[0,i,v,t] <= 1;

# (2h) Restricao de duracao da jornada
s.t. c2h{v in V, t in T diff Tprime}:
    TU * sum{i in I} x[i,0,v,t]
    + sum{i in I union {0}, j in I union {0}} x[i,j,v,t] * C[i,j]
    + sum{i in I union {0}, j in I union {0}, b in B} Sb[b] * z[i,j,b,v,t] <= TL;

# (2i) Capacidade de carga do veiculo
s.t. c2i{i in I union {0}, j in I union {0}, v in V, t in T}:
    y[i,j,v,t] <= Q * x[i,j,v,t];

# (2j) Limite inferior do balanco de carga em nos visitados
s.t. c2j_lb{j in I, v in V, t in T}:
    sum{i in I union {0}} y[j,i,v,t]
    >= sum{i in I union {0}} y[i,j,v,t] + w[j,t]
       - Q * (1 - sum{i in I union {0}} x[i,j,v,t]);

# (2m) Residuo nao excede a capacidade maxima
s.t. c2m{i in I, t in T}:
    w[i,t] <= wmax[i];

# ============================================================================
# RESTRICOES DE ACUMULO DE RESIDUOS (versao linearizada das Eqs. 2k-2l)
# ============================================================================

# (3a) Acumulo em dias intermediarios: limite inferior se nao houver coleta
s.t. waste_accum_interm{i in I, t in T: t > 1}:
    w[i,t] >= Wi[i] + w[i,t-1] - BigM * sum{j in I union {0}, v in V} x[i,j,v,t-1];

# (3a-ub1) Limite superior do acumulo (sempre)
s.t. waste_accum_interm_ub1{i in I, t in T: t > 1}:
    w[i,t] <= Wi[i] + w[i,t-1];

# (3a-ub2) Se houve coleta em t-1, residuo volta para a geracao diaria
s.t. waste_accum_interm_ub{i in I, t in T: t > 1}:
    w[i,t] <= Wi[i] + BigM * (1 - sum{j in I union {0}, v in V} x[i,j,v,t-1]);

# (3b) Acumulo no primeiro dia (ciclico): limite inferior se nao coletou no ultimo dia
s.t. waste_accum_first{i in I}:
    w[i,1] >= Wi[i] + w[i,card(T)] - BigM * sum{j in I union {0}, v in V} x[i,j,v,card(T)];

# (3b-ub1) Limite superior do acumulo para o primeiro dia (ciclico)
s.t. waste_accum_first_ub1{i in I}:
    w[i,1] <= Wi[i] + w[i,card(T)];

# (3b-ub2) Se houve coleta no ultimo dia, primeiro dia volta para geracao diaria
s.t. waste_accum_first_ub{i in I}:
    w[i,1] <= Wi[i] + BigM * (1 - sum{j in I union {0}, v in V} x[i,j,v,card(T)]);

# (3c) Residuo minimo em cada dia
s.t. waste_minimum{i in I, t in T}:
    w[i,t] >= Wi[i];

# ============================================================================
# LINEARIZACAO DO TEMPO DE SERVICO (variaveis z) - Glover (1975)
# ============================================================================

# (4a) z <= nb (se a combinacao b nao esta em j, entao z = 0)
s.t. c4a{i in I union {0}, j in I union {0}, b in B, v in V, t in T}:
    z[i,j,b,v,t] <= nb[b,j];

# (4b) z <= x (se o veiculo nao percorre i->j, entao z = 0)
s.t. c4b{i in I union {0}, j in I union {0}, b in B, v in V, t in T}:
    z[i,j,b,v,t] <= x[i,j,v,t];

# (4c) z >= nb + x - 1 (z = 1 sse nb=1 e x=1)
s.t. c4c{i in I union {0}, j in I union {0}, b in B, v in V, t in T}:
    z[i,j,b,v,t] >= nb[b,j] + x[i,j,v,t] - 1;

# ============================================================================
# CORTES DE QUEBRA DE SIMETRIA
# ============================================================================

# (6a) Ordenacao de veiculos: v so pode ser usado se v-1 for usado
s.t. c6a{v in V, t in T diff Tprime: card(V) > 1 and v > 1}:
    sum{j in I} x[0,j,v,t] <= sum{j in I} x[0,j,v-1,t];

# (6b) Veiculos iniciam com carga vazia
s.t. c6b{v in V, t in T diff Tprime}:
    sum{j in I} y[0,j,v,t] = 0;

end;
<= sum{j in I} x[0,j,v-1,t];

# (6b) Vehicles start with empty load
s.t. c6b{v in V, t in T diff Tprime}:
    sum{j in I} y[0,j,v,t] = 0;

end;<= sum{j in I} x[0,j,v-1,t];

# (6b) Vehicles start with empty load
s.t. c6b{v in V, t in T diff Tprime}:
    sum{j in I} y[0,j,v,t] = 0;

end;<= sum{j in I} x[0,j,v-1,t];

# (6b) Vehicles start with empty load
s.t. c6b{v in V, t in T diff Tprime}:
    sum{j in I} y[0,j,v,t] = 0;

end;