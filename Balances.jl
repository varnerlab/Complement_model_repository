# ----------------------------------------------------------------------------------- #
# Copyright [c] 2016 Varnerlab
# School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14853

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files [the "Software"], to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

function Balances(t,x,dxdt_vector,data_dictionary)


#-------------------------------------------
# COMPLEMENT SPECIES
#-------------------------------------------
Initiator2          = x[1];   #   Initiator complex 2, lectin pathway
C4                  = x[2];   #   C4 protein
C2                  = x[3];   #   C2 protein
C4a                 = x[4];   #   C4a protein
C4b                 = x[5];   #   C4b protein
C2a                 = x[6];   #   C2a protein
C2b                 = x[7];   #   C2b protein
C3                  = x[8];   #   C3 protein
C3b                 = x[9];   #   C3b protein
C4bC2a              = x[10];  #   Classical/Lectin C3 convertase
C3Convertase2       = x[11];  #   Alternate C3 convertase
C4bC2aC3b           = x[12];  #   Classical C5 convertase
C5Convertase2       = x[13];  #   Alternate C5 convertase
C5                  = x[14];  #   C5 protein
C5a                 = x[15];  #   C5a protein
C5b                 = x[16];  #   C5b protein
C4BP                = x[17];  #   C4 binding protein
FactorH             = x[18];  #   Factor H
C3a                 = x[19];  #   C3a protein

if (Initiator2<1e-3)
  Initiator2 = 0.0
end

if (C5<0)
  C5 = 0.0
end

if (C3<0)
  C3 = 0.0
end


#-------------------------------------------
# parameter_array
#-------------------------------------------
parameter_array = data_dictionary["PARAMETER_ARRAY"]
k_Initiator2_C4             = 1.0*parameter_array[1];
Km_Initiator2_C4            = 1.0*parameter_array[2];
k_Initiator2_C2             = 1.0*parameter_array[3];
Km_Initiator2_C2            = 1.0*parameter_array[4];

n_In2_C4                    = parameter_array[5];
n_In2_C2                    = parameter_array[6];
alpha_In2_C4                = parameter_array[7];
alpha_In2_C2                = parameter_array[8];

n_C5a                       = parameter_array[9];
n_C3a                       = parameter_array[10];
k_C4bC2a                    = parameter_array[11];
k_C3b_basal                 = parameter_array[12];
k_C3Convertase2             = parameter_array[13];
k_C5Convertase2             = parameter_array[14];
k_C5_C4bC2a_bind            = parameter_array[15];
k_C5_conversion             = parameter_array[16];
Km_C5_conversion            = parameter_array[17];
k_C5_conversion_alternate   = parameter_array[18];
Km_C5_conversion_alternate  = parameter_array[19];
k_inhibit_C4bC2a            = parameter_array[20];
k_inhibit_C3Convertase2     = parameter_array[21];
k_cat_C4bC2a                = parameter_array[22];
Km_C4bC2a                   = parameter_array[23];
k_cat_C3convertase2         = parameter_array[24];
Km_C3Convertase2            = parameter_array[25];
k_degradationC3a            = parameter_array[26];
kC5inhibit                  = parameter_array[27];
k_degradationC5a            = parameter_array[28];

#-------------------------------------------
# REACTIONS
#-------------------------------------------

# Formation Reactions
rV = zeros(14)
rV[1]  = k_Initiator2_C4*Initiator2*C4/(Km_Initiator2_C4+C4);                                             #   Lectin Pathway initiator generating C4a and C4b from C4
rV[2]  = k_Initiator2_C2*Initiator2*C2/(Km_Initiator2_C2+C2);                                             #   Lectin Pathway generating C2a and C2b from C2
rV[3]  = k_C4bC2a*C4b*C2a;                                                                                #   Formation of classical/lectin C3 convertase
rV[4]  = k_C3b_basal*C3;                                                                                  #   Formation of basal C3b and C3a from alternate pathway
rV[5]  = k_C3Convertase2*C3b;                                                                             #   Formation of alternate C3 convertase
rV[6]  = k_C5Convertase2*C3Convertase2*C3b;                                                               #   Formation of alternate C5 convertase
rV[7]  = k_C5_C4bC2a_bind*C4bC2a*C3b;                                                                     #   Formation of classical C5 convertase
rV[8] = k_C5_conversion*C4bC2aC3b*C5^n_C5a/(Km_C5_conversion^n_C5a+C5^n_C5a);                             #   Formation of C5a and C5b through classical C5 convertase
rV[9] = k_C5_conversion_alternate*C5Convertase2*C5/(Km_C5_conversion_alternate+C5);                       #   Formation of C5a and C5b through alternate C5 convertase

# Inhibition Reactions
rV[10] = k_inhibit_C4bC2a*C4b*C4BP;                                                                       #   Inhibition of C4b through C4BP
rV[11] = k_inhibit_C3Convertase2*C3Convertase2*FactorH;                                                   #   Inhibition of C3b through Factor H
rV[12] = kC5inhibit*C4bC2aC3b*C4BP;                                                                       #   Inhibition of Classical C5 Convertase by C4BP

# Formation of C3a and C3b from C3 Convertase
rV[13] = k_cat_C4bC2a*C4bC2a*C3^n_C3a/(Km_C4bC2a^n_C3a+C3^n_C3a);                                         #   Formation of C3a and C3b from Lectin C3 Convertase
rV[14] = k_cat_C3convertase2*C3Convertase2*C3/(Km_C3Convertase2+C3);                                      #   Formation of C3a and C3b from Alternate C3 Convertase

#Modified Initiator Steps -- Includes control functions
rV_modified = zeros(2)

# Setup logistics function -
control_vector = zeros(2)
# if (Initiator2 > 1e-5)
#
#   L = 1.0
#
#   alpha_C2 = alpha_In2_C2
#   alpha_C4 = alpha_In2_C4
#   k_C2 = n_In2_C2
#   k_C4 = n_In2_C4
#
#   control_vector[1] = L/(1+exp(-k_C4*(Initiator2 - alpha_C4)))
#   control_vector[2] = L/(1+exp(-k_C2*(Initiator2 - alpha_C2)))
#
# end

control_vector[1] = (Initiator2^(n_In2_C4))/(Initiator2^(n_In2_C4)+alpha_In2_C4^(n_In2_C4))
control_vector[2] = (Initiator2^(n_In2_C2))/(Initiator2^(n_In2_C2)+alpha_In2_C2^(n_In2_C2))

# @show control_vector,Initiator2
tau_C5 = 1.0
if (Initiator2>2e-5)
  tau_C5 = Initiator2;
end

rV_modified[1] = rV[1]*control_vector[1]                                                               #   Modified reaction rate for C4 conversion
rV_modified[2] = rV[2]*control_vector[2]                                                               #   Modified reaction rate for C2 conversion


# if (Initiator2 > 1e-5)
#   tau_C5 = control_vector[2]
# end

# rV_modified[1] = rV[1]*control_vector[1]
# rV_modified[1] = rV[2]*control_vector[2]

#-------------------------------------------
# MATERIAL BALANCES
#-------------------------------------------

dxdt_vector[1]   = 0;                                                             # Lectin Pathway Initiator
dxdt_vector[2]   = -rV_modified[1];                                               # C4
dxdt_vector[3]   = -rV_modified[2];                                               # C2
dxdt_vector[4]   = rV_modified[1];                                                # C4a
dxdt_vector[5]   = rV_modified[1]-rV[3];                                          # C4b
dxdt_vector[6]   = rV_modified[2]-rV[3];                                          # C2a
dxdt_vector[7]   = rV_modified[2];                                                # C2b
dxdt_vector[8]   = -rV[4]-rV[13]-rV[14];                                          # C3
dxdt_vector[9]   = rV[4]-rV[5]-rV[6]-rV[7]+rV[13]+rV[14];                         # C3b
dxdt_vector[10]  = rV[3]-rV[7]-rV[10];                                            # C4bC2a
dxdt_vector[11]  = rV[5]-rV[6]-rV[11];                                            # C3Convertase2
dxdt_vector[12]  = rV[7]-rV[12];                                                  # C4bC2aC3b
dxdt_vector[13]  = rV[6];                                                         # C5Convertase2
dxdt_vector[14]  = -1*(rV[8]+rV[9]);                                              # C5
dxdt_vector[15]  = tau_C5*(rV[8]+rV[9]-k_degradationC5a*C5a);                     # C5a
dxdt_vector[16]  = tau_C5*(rV[8]+rV[9]);                                          # C5b
dxdt_vector[17]  = -rV[10]-rV[12];                                                # C4BP
dxdt_vector[18]  = -rV[11];                                                       # FactorH
dxdt_vector[19]  = rV[4]+rV[13]+rV[14]-k_degradationC3a*C3a;                      # C3a

end
