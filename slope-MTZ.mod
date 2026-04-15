set TARGETS := {1..10}; # a set of targets
set DRONES := {"Drone1", "Drone2"}; # a set of drones
set OPERATIONS := {1..2}; # a set of operations
param n := card(TARGETS);
param subtours := 2^n - 2;
set S{i in 1..subtours} := {j in 1..n: (i div 2^(j-1)) mod 2 = 1};

#parameter
param omega_m >= 0 default 50; 
param omega_u >= 0 default 20;
param omega_p >= 0 default 1000;
param Du_max >= 0 default 7.83;
param origin_x := 0;
param origin_y := 0;
param dest_x := 4;
param dest_y := 4;
param Vm := 1;
param Vu{DRONES};
param t_charge_i{DRONES};
param A_ji {TARGETS, 1..10} >= 0;
param d_jk {j in TARGETS, k in TARGETS} = sqrt((A_ji[j, 1] - A_ji[k, 1])^2 + (A_ji[j, 2] - A_ji[k, 2])^2);
param M := 1000;

#decision variables
var x_sj {TARGETS, DRONES, OPERATIONS} binary;
var x_jf {TARGETS, DRONES, OPERATIONS} binary;
var x_j {TARGETS, DRONES, OPERATIONS} binary;
var x_jk {TARGETS, TARGETS, DRONES, OPERATIONS} binary;
var p_k {TARGETS} binary;
var h {DRONES, OPERATIONS} binary;

#contiuous variables
var d_sj {TARGETS, DRONES, OPERATIONS} >= 0;
var d_jf {TARGETS, DRONES, OPERATIONS} >= 0;
var d_LR {DRONES, OPERATIONS} >= 0;
var d_RL {DRONES, OPERATIONS} >= 0;
var d_u {DRONES, OPERATIONS} >= 0;
var d_origin >= 0;
var d_dest >= 0;
var XL_x {DRONES, OPERATIONS} >= 0;  
var XL_y {DRONES, OPERATIONS} >= 0;
var XR_x {DRONES, OPERATIONS} >= 0;  
var XR_y {DRONES, OPERATIONS} >= 0;
var subtour {TARGETS} >= 0; 
var slope >= 0;

minimize TotalOperationalCost:   
    omega_m * 
    (
        d_origin 
        + sum {o in OPERATIONS} d_LR["Drone1",o] 
        + sum {o in OPERATIONS} d_RL["Drone1",o] 
        #+ d_dest
    )
    + omega_u * (
        sum { i in DRONES, o in OPERATIONS} d_u[i,o] ) 
    + sum {k in TARGETS} omega_p * p_k[k];
    
s.t. LaunchPointY {i in DRONES, o in OPERATIONS}: 
	XL_y[i,o] = slope * XL_x[i,o];

s.t. RetrievelPointY {i in DRONES, o in OPERATIONS}: 
	XR_y[i,o] = slope * XR_x[i,o];

s.t. ValidLaunchPoint {i in DRONES, o in OPERATIONS}:
    sum {j in TARGETS} x_sj[j, i, o] <= 1;
    
s.t. TargetVisit {j in TARGETS, i in DRONES, o in OPERATIONS}:
    x_j[j, i, o] = x_sj[j, i, o] + sum {k in TARGETS: k != j} x_jk[k,j, i, o];

s.t. DroneMovement {j in TARGETS, i in DRONES, o in OPERATIONS}:
    x_j[j, i, o] = x_jf[j, i, o] + sum {k in TARGETS: k != j} x_jk[j,k, i, o]; 

s.t.  ValidRetrievalPoint {i in DRONES, o in OPERATIONS}:
    sum {j in TARGETS} x_jf[j, i, o] <= 1;

s.t. DronesTotalDistance {i in DRONES, o in OPERATIONS}:
    d_u[i, o] = 
        sum {j in TARGETS} x_sj[j, i, o] * d_sj[j, i, o] + 
        sum {j in TARGETS, k in TARGETS: j != k} x_jk[j,k, i, o] * d_jk[j, k] +
        sum {j in TARGETS} x_jf[j, i, o] * d_jf[j, i, o];
              
s.t.  DronesBatteryCapacity {i in DRONES, o in OPERATIONS}:
    d_u[i, o] <= Du_max;
    
s.t. dLR_def {i in DRONES, o in OPERATIONS}:
	d_LR[i,o] >= (d_u[i, o] / Vu[i]) * Vm;  #ori =
	
s.t. LR_CoordGap {i in DRONES, o in OPERATIONS}:
    XR_x[i,o] - XL_x[i,o] >= d_LR[i,o] / sqrt(1 + slope^2);

s.t. RL_CoordGap {i in DRONES, o in OPERATIONS: o > 1}:
    XL_x[i,o] - XR_x[i,o-1] >= d_RL[i,o] / sqrt(1 + slope^2);

s.t. RechargeGap {i in DRONES, o in OPERATIONS: o > 1}:
    d_RL[i,o] >= Vm * t_charge_i[i];


s.t. SubTourElimination {j in TARGETS, k in TARGETS, i in DRONES, o in OPERATIONS : j != k}:
    subtour[j] - subtour[k] + M * x_jk[j,k, i, o] <= M - 1;

s.t. VisitAllTargets {k in TARGETS}:
    sum {i in DRONES, o in OPERATIONS} 
        (x_sj[k, i, o] + sum {j in TARGETS: j != k} x_jk[j, k, i, o]) = 1 - p_k[k];
        
s.t. DIST1 {j in TARGETS, i in DRONES, o in OPERATIONS}: 
    d_sj[j, i, o] = sqrt((XL_x[i, o] - A_ji[j, 1])^2 + (XL_y[i, o] - A_ji[j, 2])^2);

s.t.  DIST3 {j in TARGETS, i in DRONES, o in OPERATIONS}:
    d_jf[j, i, o] =  sqrt((A_ji[j, 1] - XR_x[i, o])^2 + (A_ji[j, 2] - XR_y[i, o])^2);
      
s.t. DIST6:
    d_origin = sqrt((origin_x - XL_x["Drone1", 1])^2 + (origin_y - XL_y["Drone1", 1])^2);


s.t. UsedOpFromLaunch {i in DRONES, o in OPERATIONS}:
    sum {j in TARGETS} x_sj[j,i,o] = h[i,o];

s.t. UsedOpFromRetrieve {i in DRONES, o in OPERATIONS}:
    sum {j in TARGETS} x_jf[j,i,o] = h[i,o];

s.t. OperationOrderSymmetry {i in DRONES, o in OPERATIONS: o < card(OPERATIONS)}:
    h[i,o] >= h[i,o+1];
    


    