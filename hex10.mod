# SETS & PARAMETERS
set Q; set R;
set GRID within {Q,R};
set TARGETS within GRID;
set DRONES := {"Drone1","Drone2"};
set OPERATIONS := {1..2};

param TM integer >= 2;     
param TDmax integer >= 1;   
set TMSET := 1..TM;
set TDSET := 1..TDmax;
set M_ARC := {m in TDSET: m < TDmax};
  
param d_step >= 0 default 1;             
param Vu{DRONES} >= 0;                
param Vm := 1; # mothership speed (distance per unit of time)                         
param omega_m >= 0 default 50;           
param omega_u >= 0 default 20;  
param omega_p >= 0 default 1000;           
param t_charge{DRONES} integer >= 0; 
param Delta_m := d_step / Vm; # Duration of 1 TM step (time per step)
param Delta_u{i in DRONES} := d_step / Vu[i];  # Duration of 1 TD step (time per step)
param charge_steps_TM{i in DRONES} := ceil( t_charge[i] / Delta_m );  # TM step
param t_max{DRONES} integer >= 0; 
param origin_q in Q;  
param origin_r in R;
param dest_q in Q;  
param dest_r in R;
param M := 1000;

# Neighbors 
set NEIGHBORS {(q1,r1) in GRID} within GRID :=
   setof{(q2,r2) in GRID :
    ( q1 mod 2 = 1 and
      (
        (q2 = q1-1 and r2 = r1+1) or
        (q2 = q1   and r2 = r1+1) or
        (q2 = q1+1 and r2 = r1+1) or
        (q2 = q1+1 and r2 = r1  ) or
        (q2 = q1   and r2 = r1-1) or
        (q2 = q1-1 and r2 = r1  )
      )
    )
    or
    ( q1 mod 2 = 0 and
      (
        (q2 = q1-1 and r2 = r1  ) or
        (q2 = q1   and r2 = r1+1) or
        (q2 = q1+1 and r2 = r1  ) or
        (q2 = q1+1 and r2 = r1-1) or
        (q2 = q1   and r2 = r1-1) or
        (q2 = q1-1 and r2 = r1-1)
      )
    )
  } (q2,r2);

#Mothership ARCS
set ARCS within {GRID, GRID} := 
  setof {(q1,r1) in GRID, (q2,r2) in NEIGHBORS[q1,r1]} (q1,r1,q2,r2);
set DEST_STAY := {(dest_q, dest_r, dest_q, dest_r)};
set ARCS_STAY := ARCS union DEST_STAY;
 
# VARIABEL (Mothership)
var ypos {GRID, t in TMSET} binary; 
var yarc {(q1,r1,q2,r2) in ARCS_STAY, t in TMSET diff {TM} } binary;
var yvisit {(q,r) in TARGETS} binary; # the mothership visits target   
var m_arrive_t {t in TMSET} binary;   # 1 if first time at (dest_q,dest_r) is t

 
# VARIABLE - Operation Synchronization (at mothership time)
var s {i in DRONES, o in OPERATIONS, (q,r) in GRID, t in TMSET} binary;  # start op o
var f {i in DRONES, o in OPERATIONS, (q,r) in GRID, t in TMSET} binary;  # finish op o
var t_l {i in DRONES, o in OPERATIONS} integer >= 0; 
var t_r {i in DRONES, o in OPERATIONS} integer >= 0; 

# VARIABEL - Drone
var z {i in DRONES, o in OPERATIONS, GRID, m in TDSET} binary; # drone occupancy on m
var x {DRONES, OPERATIONS, TARGETS} binary;
var xarc {i in DRONES, o in OPERATIONS, ARCS, m in TDSET diff {TDmax}} binary; 
#var xarc{i in DRONES, o in OPERATIONS, (q1,r1,q2,r2) in ARCS, m in TDSET diff {TDmax}} binary;
var t_end {i in DRONES, o in OPERATIONS} integer >= 0 <= TDmax; # panjang langkah drone (akhir operasi, dalam m)
var x_end {i in DRONES, o in OPERATIONS, m in TDSET} binary;   # pilih 1 m sbg akhir

# active/inactive operation
var h {i in DRONES, o in OPERATIONS} binary;
var p { (q,r) in TARGETS } binary;

# objective function
minimize Total_Cost:
    (omega_m * d_step *
      sum{ (q1,r1,q2,r2) in ARCS, t in TMSET: t < TM } yarc[q1,r1,q2,r2,t]
  + omega_u * d_step *
      sum{ i in DRONES, o in OPERATIONS, (q1,r1,q2,r2) in ARCS, m in TDSET: m < TDmax }
          xarc[i,o,q1,r1,q2,r2,m]
  + omega_p * sum{ (q,r) in TARGETS } p[q,r])*1.12;

s.t. M_Start:
    ypos[origin_q, origin_r, 1] = 1;

s.t. M_End:
    ypos[dest_q, dest_r, TM] = 1;

s.t. M_FlowOut { (q,r) in GRID, t in TMSET: t < TM }:
  ypos[q,r,t] = sum{ (q2,r2) in GRID: (q,r,q2,r2) in ARCS_STAY } yarc[q,r,q2,r2,t];


s.t. M_FlowIn { (q,r) in GRID, t in TMSET: t < TM }:
  ypos[q,r,t+1] = sum{ (q1,r1) in GRID: (q1,r1,q,r) in ARCS_STAY } yarc[q1,r1,q,r,t];

s.t. M_OnePosPerT { t in TMSET }:
    sum{(q,r) in GRID} ypos[q,r,t] = 1; 
    
s.t. StartAtMothership {i in DRONES, o in OPERATIONS, (q,r) in GRID, t in TMSET}:
    s[i,o,q,r,t] <= ypos[q,r,t];

s.t. FinishAtMothership {i in DRONES, o in OPERATIONS, (q,r) in GRID, t in TMSET}:
    f[i,o,q,r,t] <= ypos[q,r,t];

s.t. Def_tl {i in DRONES, o in OPERATIONS}:
    t_l[i,o] = sum{(q,r) in GRID, t in TMSET} t * s[i,o,q,r,t];

s.t. Def_tr {i in DRONES, o in OPERATIONS}:
    t_r[i,o] = sum{(q,r) in GRID, t in TMSET} t * f[i,o,q,r,t];

s.t. Valid_Launch_Point {i in DRONES, o in OPERATIONS}:
    sum { (q,r) in GRID, t in TMSET } s[i,o,q,r,t] = h[i,o];

s.t. Valid_Retrieve_Point {i in DRONES, o in OPERATIONS}:
    sum { (q,r) in GRID, t in TMSET } f[i,o,q,r,t] = h[i,o];  
    
s.t. D_EndExactlyOneCell_LB {i in DRONES, o in OPERATIONS, m in TDSET}:
    sum{(q,r) in GRID} z[i,o,q,r,m] >= x_end[i,o,m];

s.t. D_AtMostOneCellPerM {i in DRONES, o in OPERATIONS, m in TDSET}:
    sum{(q,r) in GRID} z[i,o,q,r,m] <= 1;

s.t. D_Propagate {i in DRONES, o in OPERATIONS, (q,r) in GRID, m in TDSET: m < TDmax}:
   z[i,o,q,r,m+1] = sum{(q1,r1) in NEIGHBORS[q,r]} xarc[i,o,q1,r1,q,r,m];

s.t. NoArcAtEnd {i in DRONES, o in OPERATIONS, (q1,r1,q2,r2) in ARCS, m in TDSET: m < TDmax}:
  xarc[i,o,q1,r1,q2,r2,m] <= 1 - x_end[i,o,m];


s.t. OutflowUB {i in DRONES, o in OPERATIONS, (q,r) in GRID, m in TDSET: m < TDmax}:
  sum{(q2,r2) in NEIGHBORS[q,r]} xarc[i,o,q,r,q2,r2,m] <= z[i,o,q,r,m];

s.t. OutflowLB {i in DRONES, o in OPERATIONS, (q,r) in GRID, m in TDSET: m < TDmax}:
  sum{(q2,r2) in NEIGHBORS[q,r]} xarc[i,o,q,r,q2,r2,m] >= z[i,o,q,r,m]-x_end[i,o,m];


s.t. D_StartOcc {i in DRONES, o in OPERATIONS, (q,r) in GRID}:
    z[i,o,q,r,1] = sum{t in TMSET} s[i,o,q,r,t];

s.t. EndOne {i in DRONES, o in OPERATIONS}:
    sum{m in TDSET} x_end[i,o,m] = h[i,o];
    
s.t. MendDef {i in DRONES, o in OPERATIONS}:
    t_end[i,o] = sum{m in TDSET} m * x_end[i,o,m]; 
 
s.t. FinishLink_LB {i in DRONES, o in OPERATIONS, (q,r) in GRID, m in TDSET}:
    sum{t in TMSET} f[i,o,q,r,t] >= z[i,o,q,r,m] - (1 - x_end[i,o,m]);


s.t. Link_Target_To_Occupancy {i in DRONES, o in OPERATIONS, (q,r) in TARGETS}:
    x[i,o,q,r] = sum{m in TDSET} z[i,o,q,r,m]; #try =
      
s.t. CoverAll_Soft{(q,r) in TARGETS}:
    yvisit[q,r] + sum{i in DRONES, o in OPERATIONS} x[i,o,q,r] >= 1 - p[q,r];

s.t. Yvisit_Upper { (q,r) in TARGETS }:
    yvisit[q,r] <= sum { t in TMSET } ypos[q,r,t];

s.t. TimeSync_Physical {i in DRONES, o in OPERATIONS}:
    ( t_r[i,o] - t_l[i,o] ) * Delta_m  >=  t_end[i,o] * Delta_u[i];

s.t. NextOp_Launch_TM {i in DRONES, o in OPERATIONS: o < card(OPERATIONS)}:
  t_l[i,o+1] >= t_r[i,o] + charge_steps_TM[i];
  
s.t. ArriveFlag1 {t in TMSET}: m_arrive_t[t] <= ypos[dest_q, dest_r, t]; 
 

s.t. ArriveUnique: sum {t in TMSET} m_arrive_t[t] = 1;

s.t. RetrieveBeforeArrival {i in DRONES, o in OPERATIONS}:
    t_r[i,o] <= sum {t in TMSET} t * m_arrive_t[t];

s.t. LaunchBeforeArrival {i in DRONES, o in OPERATIONS}:
    t_l[i,o] <= sum {t in TMSET} (t-1) * m_arrive_t[t];  
   
s.t. ArrivePrefix {t in TMSET}:
    sum {tau in TMSET: tau <= t} m_arrive_t[tau] >= ypos[dest_q, dest_r, t]; 

     