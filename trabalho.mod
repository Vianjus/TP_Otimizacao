# ============================================================================
# Periodic Vehicle Routing and Capacitated Facility Location Problem
# for Municipal Solid Waste Collection - MILP Exact Method
# ============================================================================

set I;                          # Collection points (excluding depot)
set B;                          # Bin combinations  
set V;                          # Vehicles
set T;                          # Days in planning horizon

param Q;                        # Vehicle capacity (m3)
param C{i in I union {0}, j in I union {0}};  # Travel time between points (min)
param Sb{b in B};              # Service time of bin combination b (min)
param TU;                       # Unloading time at depot (min)
param Wi{i in I};              # Daily waste generation at point i (m3/day)
param CAPb{b in B};            # Capacity of bin combination b (m3)
param CINb{b in B};            # Installation & maintenance cost of bin combination b (US$)
param CCV;                      # Cost of vehicle operations (US$/min)
param TL;                       # Working day length for drivers (min)
param BigM;                     # Big M for linearization

# ============================================================================
# DECISION VARIABLES
# ============================================================================

# Routing variables
var x{i in I union {0}, j in I union {0}, v in V, t in T}, binary;
var y{i in I union {0}, j in I union {0}, v in V, t in T}, >= 0;

# Waste accumulation variables
var w{i in I, t in T}, >= 0;
var wmax{i in I}, >= 0;

# Bin allocation variables
var nb{b in B, i in I union {0}}, binary;

# Linearization variables for service time
var z{i in I, b in B, v in V, t in T}, >= 0;

# ============================================================================
# OBJECTIVE FUNCTION
# ============================================================================

minimize Total_Cost:
    sum{b in B, i in I} CINb[b] * nb[b,i]
    + CCV * sum{t in T, v in V} (
        TU * sum{i in I} x[i,0,v,t]
        + sum{i in I union {0}, j in I union {0}} x[i,j,v,t] * C[i,j]
        + sum{i in I, b in B} Sb[b] * z[i,b,v,t]
      );

# ============================================================================
# CONSTRAINTS
# ============================================================================

# (2a) No bin combination at depot
s.t. c2a{b in B}:
    nb[b,0] = 0;

# (2b) Exactly one bin combination per collection point
s.t. c2b{i in I}:
    sum{b in B} nb[b,i] = 1;

# (2c) Bin capacity >= maximum waste accumulated
s.t. c2c{i in I}:
    sum{b in B} CAPb[b] * nb[b,i] >= wmax[i];

# (2d) No self-loops
s.t. c2d{i in I union {0}, v in V, t in T}:
    x[i,i,v,t] = 0;

# (2f) Flow conservation
s.t. c2f{j in I union {0}, v in V, t in T}:
    sum{i in I union {0}} x[i,j,v,t] = sum{i in I union {0}} x[j,i,v,t];

# (2g) At most one route per vehicle per day
s.t. c2g{v in V, t in T}:
    sum{i in I} x[0,i,v,t] <= 1;

# (2i) Vehicle load capacity
s.t. c2i{i in I union {0}, j in I union {0}, v in V, t in T}:
    y[i,j,v,t] <= Q * x[i,j,v,t];

# (2j) Flow balance for load (subtour elimination + load tracking)
s.t. c2j{j in I, v in V, t in T}:
    sum{i in I union {0}} y[i,j,v,t] + w[j,t] 
    <= sum{i in I union {0}} y[j,i,v,t] 
       + Q * (1 - sum{i in I union {0}} x[i,j,v,t]);

# (2m) Waste does not exceed bin capacity
s.t. c2m{i in I, t in T}:
    w[i,t] <= wmax[i];

# Basic waste accumulation: each point accumulates waste daily
s.t. basic_waste{i in I, t in T}:
    w[i,t] >= Wi[i];

# ============================================================================
# LINEARIZATION FOR SERVICE TIME (z variables)
# ============================================================================

# z[i,b,v,t] = 1 iff point i is visited by vehicle v on day t AND bin b is installed at i
# (4a) z <= nb (if bin b not at i, then z = 0)
s.t. c4a{i in I, b in B, v in V, t in T}:
    z[i,b,v,t] <= nb[b,i];

# (4b) z <= x (if point i not visited, then z = 0)  
s.t. c4b{i in I, b in B, v in V, t in T}:
    z[i,b,v,t] <= sum{j in I union {0}} x[j,i,v,t];

# (4c) z >= nb + x - 1 (z = 1 iff both nb = 1 and visited)
s.t. c4c{i in I, b in B, v in V, t in T}:
    z[i,b,v,t] >= nb[b,i] + sum{j in I union {0}} x[j,i,v,t] - 1;

# ============================================================================
# SYMMETRY BREAKING CUTS
# ============================================================================

# (6a) Vehicle ordering: vehicle v can only be used if vehicle v-1 is used
s.t. c6a{v in V, t in T : v > 1}:
    sum{j in I} x[0,j,v,t] <= sum{j in I} x[0,j,v-1,t];

# (6b) Vehicles start with empty load
s.t. c6b{v in V, t in T}:
    sum{j in I} y[0,j,v,t] = 0;

end;
