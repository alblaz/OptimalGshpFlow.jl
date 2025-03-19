### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bc6d9b90-734e-4d9d-be0f-104a3b99aea2
using PlutoUI

# ╔═╡ 8857ee4a-9620-4071-ae78-6e63ea0b7f0e
using CoolProp

# ╔═╡ a9eddcc7-ac94-483e-8622-e6c9dca8351b
using Unitful

# ╔═╡ 0405ca05-0ca9-4b6c-ac93-91e5f50cab91
using Markdown

# ╔═╡ c691da35-ff83-4dc6-af9f-2be03d95187a
using InteractiveUtils

# ╔═╡ f6b9843f-bb5e-4275-bea2-05e9d594207f
begin
    using WGLMakie
    set_theme!(theme_dark())
end

# ╔═╡ 1c61fab6-f56b-460b-870f-8e3a83337312
using ModelingToolkit, NonlinearSolve

# ╔═╡ fbf3c417-8b31-4256-8981-84b0f0ce34b0
using LinearAlgebra

# ╔═╡ fa00c042-0b93-4635-ac75-4513e61671ec
using Nemo

# ╔═╡ 8b583195-b96d-4cd6-8bb5-f09512c38ea6
md""" # TDL:
- Connect the inlets and outlets of the 2 HX (BH/Evaporator)
- Use the outlet coindensor enthalpy as an input to the evaporator
- Have dynamic HTCs -> reuse the previous wok for the BH & Use Alberto's correlations for the evaporator
- Solve with the HTCs
- Implementing other types of collector
- etc.
"""

# ╔═╡ 4b6d717a-edb6-4292-bcb7-5dc44e560637
md"## Borehole Heat Exchanger"

# ╔═╡ 7e0fdfd5-dc29-4cb1-9ee5-82442ec8933c
md"""
Length: $(@bind Lbh_ Scrubbable(10.0:5:1000.0, default = 300)) m"""

# ╔═╡ 81a7c7ba-9cba-4b50-98c8-e14f5edf64d3
md"""
Borehole diameter: $(@bind Dbh_ Scrubbable(5:0.5:30.0, default = 11.5)) cm"""

# ╔═╡ ac5f77da-a6b7-4be1-8fb6-523c10f2a147
md"""
Type of collector: $(@bind collector Select(["U pipe" => "U pipe","Double U pipe" => "Double U pipe", "Coaxial" => "Coaxial"]))"""

# ╔═╡ 9e09fbba-ecdf-4615-8baf-f25a7cafb195
md"""
Pipe diameter: $(@bind Dp_ Scrubbable(5:0.5:30.0, default = 11.5)) cm"""

# ╔═╡ b754c0a6-4579-40f2-b03c-31142ef5bb2f
md"""
UA pipe to wall: $(@bind Upw_ Scrubbable(1.0:0.1:10.0, default = 5.0)) W/m2/K"""

# ╔═╡ 2dffa4af-cb4e-42f0-be0e-417734d119d4
md"""
UA pipe to pipe: $(@bind Upp_ Scrubbable(1.0:0.1:10.0, default = 5.0)) W/m2/K"""

# ╔═╡ 0a284427-2f1e-4d6f-b3d9-32b0b54c4564
md" ## Secondary fluid"

# ╔═╡ 7182187c-c0a3-4b8b-bd7e-06718192666a
md"""
Secondary fluid: $(@bind fluid Select(["Water" => "Water", "MEA" => "Ethanol - Water", "MEG" => "Ethylene Glycol - Water"]))"""

# ╔═╡ 4af116ad-3eac-479d-8c44-2e90ab2021c2
md"""
Temperature difference fluid-superheat: $(@bind ΔT_brine_ Scrubbable(1.0:0.5:15.0; default=5.0)) °C
"""

# ╔═╡ 89c842e7-73ac-4258-83e4-a98ca23e5a12
md"## Evaporator"

# ╔═╡ 937426c1-8075-45a7-9d92-8a729b1ac251
md"""
Length: $(@bind L_ Scrubbable(0.0:0.05:10.0, default = 1.0)) m"""

# ╔═╡ 436b93eb-fc29-4529-9e6d-2aa16e808f1c
md"""
Perimeter: $(@bind P_ Scrubbable(1.0:0.5:100.0, default = 10.0)) m"""

# ╔═╡ 64756cb4-776f-4512-8f2f-fc11b9715859
md"""
U evaporation + forced convection fluid: $(@bind U_A_ Scrubbable(500:5:1500.0, default = 500.0)) W/m2/K"""

# ╔═╡ 8e6ffc58-6d02-4707-8504-4c917281e29f
md"""
U refrigerant gas + forced convection fluid: $(@bind U_B_ Scrubbable(200:1:500.0, default = 200.0)) W/m2/K"""

# ╔═╡ 30d089ba-e823-408c-9073-50b3f37a1c4c
md"## Operation"

# ╔═╡ 4937015c-8461-4cfe-8c83-51c3544f3abf
md"""
Refrigerant mass flow: $(@bind ṁ_r_ Scrubbable(0.0:0.05:10.0, default = 1.0)) kg/s"""

# ╔═╡ 25b4a88c-5c58-4d67-8d0f-ecbad5f5afae
md"## Refrigerant"

# ╔═╡ 3fc0ed88-9d4d-45c3-a0f2-48e38e9caa43
md"""
Refrigerant: $(@bind refrigerant Select(["R404A" => "R404A", "R134A" => "R134A"]))"""

# ╔═╡ 3703ecf4-8a7d-44f7-a59a-f429ad768205
md" ### Evaporator"

# ╔═╡ 30e8753a-b57c-4382-910d-fcaf08bbe6b3
md"""
Superheating degree: $(@bind ΔT_SH_ Scrubbable(1.0:0.1:10.0; default=7.0)) K
"""

# ╔═╡ 4bf19bc9-29f9-493f-830c-a6beb30fc8d6
md" ### Condensor"

# ╔═╡ ed6c8661-293b-4a95-a714-1fcf498883fa
md"""
Subcooling degree: $(@bind ΔT_SC_ Scrubbable(1.0:0.1:10.0; default=5.0)) K
"""

# ╔═╡ 1b08d63a-db57-44c7-b1ec-10acbcaa7951
md""" # Solving the system"""

# ╔═╡ d4aba512-dd8c-49a5-a7b1-eccf7ca24eab
@variables ṁ1 ṁ2 

# ╔═╡ 9902671a-a681-4b14-9e07-424eb673eef1
@parameters cp1 cp2 

# ╔═╡ d5cbb87a-d35d-4aff-b213-c89ed193619e
@variables LA LB x

# ╔═╡ c7bd0055-8fbe-4861-9284-1e08d89bcad3
@parameters U_A U_B

# ╔═╡ 1f23f14b-69e4-4366-acd8-8ded3a6ca847
@parameters P L

# ╔═╡ b6d31712-03e2-4981-95c4-d3621f9d0ed6
@parameters Tsat hdew h1iA T2iA T1oB T2iB

# ╔═╡ 7dd35a9e-12d1-4728-bae4-49c909cd7902
C1 = ṁ1*cp1

# ╔═╡ 8781d78b-1308-496f-8e80-aa1cd334dd12
C2 = ṁ2*cp2

# ╔═╡ 39ca3785-1a2f-4f1f-8877-76c9841ce098
md" ## This includes the assumption T2(LA) = T2i."

# ╔═╡ 7274b23c-f22e-4b2a-86b1-5976d25c363c
h1A = h1iA - C2/ṁ1*(Tsat - T2iA)*exp(-U_A*P*LA/C2)*exp((U_A*P*x/C2)-1)

# ╔═╡ b146988c-40cb-4576-97f1-32950538af16
T2A = Tsat - (Tsat - T2iA) * exp(-U_A*(LA-x)*P/C2)

# ╔═╡ ca92e844-5285-47de-a83d-5ae5c872e44b
T1A = Tsat

# ╔═╡ 2a4fca86-f5f9-4e0d-8a89-8898ea190cca
md" ### Sanity checks"

# ╔═╡ af00af29-d180-425b-88d5-92a9c8d86b79
substitute(T1A, x => 0.0)

# ╔═╡ 08d860a7-c75a-4fed-a752-63894211aee6
substitute(T1A, x => LA)

# ╔═╡ be4aa932-54b7-4017-9c06-59104b11aea0
substitute(T2A, x => 0.0)

# ╔═╡ cc07e1f9-8877-43e8-99c3-2d293bc00427
substitute(T2A, x => LA)

# ╔═╡ 1481a416-efa8-4ef9-9866-c159ce0f3e75
md""" ## This includes the assumptions:
- T1(LB) = T1iB.
- T2(LB) = T2oB."""

# ╔═╡ 92665c1f-eb3a-4c13-9f8a-5eed60e54835
T1B = T1oB + (T1oB-T2iB) * C2/(C1+C2)*(exp(P*U_B*(1/C1+1/C2)*(LB-x))-1)

# ╔═╡ ee30d6f9-b902-4d87-9ce7-f245573642ee
T2B = T2iB + (T1oB-T2iB) * C1/(C1+C2)*(exp(P*U_B*(1/C1+1/C2)*(LB-x))-1)

# ╔═╡ c87f7981-ccf7-4b22-ada3-7e7e944f1bf4
md" ### Sanity checks"

# ╔═╡ 744924ad-810b-4bd0-9417-165f6a48ed9a
substitute(T1B, x => 0.0)

# ╔═╡ e89934b8-49a3-4e32-8bb1-fdc6021163fa
substitute(T1B, x => LB)

# ╔═╡ 007afd7d-7941-4bda-99a4-fe02d6b1d0be
substitute(T2B, x => 0.0)

# ╔═╡ 3f1d4587-27fe-458a-8765-7465d29d6688
substitute(T2B, x => LB)

# ╔═╡ 470dccc9-6eca-4275-a532-d03889791b06
md"""
# Boundary conditiions

The following BC are already baked into the previous profiles:
- h1iA, the inlet enthalpy in the section
- T1oB, the outlet temperature from the superheated zone of the evaporator
- T2iB, the inlet temperature of the brine in the superheated zone of the evaporator
"""

# ╔═╡ 79558ba6-4469-4a40-966b-490f6d8e20f0
eq_cont_1 = hdew ~ substitute(h1A, x => LA)

# ╔═╡ 78d75ba8-0c50-48ff-8c4c-1f75498eb109
eq_cont_2 = substitute(T2A, x => LA) ~ substitute(T2B, x => 0)

# ╔═╡ 14ee0cbd-81e0-4da0-8b81-8020d8c0175e
eq_sat_T1 = Tsat ~ substitute(T1B, x => 0)

# ╔═╡ 5b8e04d3-0080-4755-9f48-526871e0e434
eq_length = L ~ LA + LB

# ╔═╡ 84adc973-3aad-4ed1-aab3-4ff3142da7d8
@named two_zone = NonlinearSystem([eq_cont_2, eq_sat_T1, eq_cont_1, eq_length], [LA, LB, ṁ2, T2iA], [L, Tsat, U_A, ṁ1, U_B, h1iA, hdew, T2iB, T1oB, cp1, cp2, P])


# ╔═╡ e8bcc3d7-fc8f-4960-b63b-90e769502f0e
tz = simplify(two_zone) |> complete

# ╔═╡ 4656a6dc-e314-4110-8d8c-e20ece6417bd
Tsat - (Tsat - T2iA) * exp(-U_A*(LA-x)*P/C2)

# ╔═╡ b0f357cc-f0ad-4a28-912f-9ee0517d49bc
T2A_n = build_function(T2A, x, LA, T2iA, Tsat, P, cp2, U_A, ṁ2, expression = Val{false})

# ╔═╡ 84ed2576-78fc-4779-a964-85bfbde1e01c
T1oB + (T1oB-T2iB) * C2/(C1+C2)*(exp(P*U_B*(1/C1+1/C2)*(LB-x))-1)

# ╔═╡ 39120525-093c-494e-8494-da3bd900b650
T1B_n = build_function(T1B, x, LB, T1oB, T2iB, P, cp1, cp2, U_B, ṁ1, ṁ2, expression = Val{false})

# ╔═╡ ba97eec7-b2fe-4c9c-8b8d-c6f20cde17b3
T2iB - (T1oB-T2iB) * C2/(C1+C2)*(exp(P*U_B*(1/C1+1/C2)*(LB-x))-1)

# ╔═╡ 4fe5b696-dc7f-4d52-b4e7-afbe9e9d7ca6
T2B_n = build_function(T2B, x, LB, T1oB, T2iB, P, cp1, cp2, U_B, ṁ1, ṁ2, expression = Val{false})

# ╔═╡ d77ec1b5-ebbb-463d-ae71-9cd880e5dd1a
md"## Borehole"

# ╔═╡ c6acaf2a-5f42-4932-acb1-632ad13689e5
@variables C1bh C2bh T1_in Ts Lbh λ1 λ2 γ1 γ2

# We assume that our homogeneous (shifted) solution is:
#   T̃₁(x) = C1*exp(λ1*x) + C2*exp(λ2*x)
#   T̃₂(x) = C1*γ1*exp(λ1*x) + C2*γ2*exp(λ2*x)
#
# The boundary conditions we impose are:
# (i) At x=0:   T1(0)=T1_in  →  T̃₁(0)=T1_in - Ts, i.e.,  C1 + C2 = T1_in - Ts.
# (ii) At x=L:  Impose T1(L) = T2(L), i.e.,
#        T̃₁(L) = T̃₂(L)
#        which gives: C1*exp(λ1*L) + C2*exp(λ2*L) = C1*γ1*exp(λ1*L) + C2*γ2*exp(λ2*L).

# Write the equations:

# ╔═╡ cea2d602-810a-4dc4-a0fe-0fb76d680b58
eq1 = C1bh + C2bh ~ T1_in - Ts

# ╔═╡ 0610a6d3-7b33-415c-bf55-9c16750785e7
eq2 = C1bh*(γ1 - 1)*exp(λ1*Lbh) + C2bh*(γ2 - 1)*exp(λ2*Lbh) ~ 0

# Solve the system for C1 and C2.

# ╔═╡ e8887d73-60e3-430a-b998-e864e7e4bcf3
sol_C = symbolic_linear_solve([eq1, eq2], [C1bh, C2bh])

# ╔═╡ 5bf5942e-46b4-4716-8e57-062a54d436e0
sol_C

# ╔═╡ 3f01b035-8fb7-4031-8a11-ef745c83606b
@variables UA m1 cp1bh Ue1 UA2 m2 cp2bh Ue2  # UA, etc. for fluids (note: UA appears for both)
# In our formulation we assume that:
#   P1 = UA/(m1*cp1)
#   Pe1 = Ue1*UA/(m1*cp1)
#   P2 = UA/(m2*cp2)
#   Pe2 = Ue2*UA/(m2*cp2)

# ╔═╡ 56e7aae4-89ae-4bb9-90a4-1ebbd6300d12
P1   = UA/(m1*cp1)

# ╔═╡ 97fd0f42-4e1e-453f-a239-925003eb103c
Pe1  = Ue1*UA/(m1*cp1)

# ╔═╡ 1c648223-2e60-4d35-83a8-d5488f79884b
P2   = UA/(m2*cp2)

# ╔═╡ 12b47c63-28d1-4a78-a733-23e10dbfe214
Pe2  = Ue2*UA/(m2*cp2)

# Next, the coefficient matrix for the homogeneous system (in shifted variables) is:
#   A = [ - (P1 + Pe1)     P1
#         - P2         P2 + Pe2 ]

# ╔═╡ 690e2521-eda6-4037-a7c3-69c48f82a25a
@variables λ  # eigenvalue placeholder

# ╔═╡ 6b85b101-74e2-4d27-9c27-1fff81802c53
A_mat = [ -(P1 + Pe1)   P1;
          -P2    P2 + Pe2 ]

# The characteristic equation is det(A - λI)=0:

# ╔═╡ 15752913-4b37-4de4-9337-86ace6dcf64f
I2 = Num.([1 0;
          0 1])

# The characteristic equation is det(A - λI)=0:

#I2 = Matrix{Num}(I, 2, 2)   # 2x2 identity matrix (Symbolics-compatible)

# ╔═╡ 3d776a5a-a56f-4094-a15e-e2e275f8ffe5
# Solve for λ symbolically:

# ╔═╡ 72c892e4-64fe-44b5-bcfb-0673289be93f
sol_λ = Symbolics.symbolic_solve(det(A_mat - λ*I2) ~ 0, λ)

# ╔═╡ 20b6adea-53c7-4d6c-91d4-9784fd6d8113
sλ1 = sol_λ[1]

# ╔═╡ ef63e77a-f931-4704-9774-dee7f37ee9eb
sλ2 = sol_λ[2]

# ╔═╡ 56016b33-a3c6-4f3b-a6e9-37687349d984
sγ1 = (P1 + Pe1 + sλ1) / P1

# ╔═╡ 39f5e272-fe0b-456b-a3f3-d268d1f2756a
sγ2 = (P1 + Pe1 + sλ2) / P1

# ╔═╡ 40a43ccc-1fad-440a-b033-66061b4dffa9
sol_C_full = [substitute(s, Dict(λ1 => simplify(sλ1),
                                        λ2 => simplify(sλ2),
                                        γ1 => simplify(sγ1),
                                        γ2 => simplify(sγ2))) for s in sol_C]

# ╔═╡ 1bc06242-a9e0-4836-a061-4aa896eb11b7
simplify.(sol_C_full)

# ╔═╡ 8d53b90c-124f-433a-a9d3-87c83fddfa22
sC1 = sol_C_full[1]

# ╔═╡ 2b9a45fb-05b1-4b54-9c96-88289d79aee8
sC2 = sol_C_full[2]

# ╔═╡ 01ad16df-5b68-40a0-9e5e-16fea88ebc8e
@variables xbh

# ╔═╡ 77caced1-75ee-42ba-8483-59cf3d10fa14
sT̃1 = sC1*exp(sλ1*xbh) + sC2*exp(sλ2*xbh)

# ╔═╡ 043640e4-5b45-4e91-a0e1-c286be1b17c1
sT1 = sT̃1 + Ts

# ╔═╡ 1f43adb4-7812-45e6-b4fb-a89f207cceb9
sT̃2 = sC1*sγ1*exp(sλ1*xbh) + sC2*sγ2*exp(sλ2*xbh)

# ╔═╡ 840b2e1c-922a-4664-a622-05a4da45a3ca
sT2 = sT̃2 + Ts

# ╔═╡ 3b34cf1d-3390-4ad5-9e69-b56dca416dbd
α = 1-(substitute(sT̃2, xbh => 0.0) / substitute(sT̃1, xbh => 0.0)) |> simplify

# ╔═╡ b7f99ca5-5531-4bfb-ab80-0b8513578d63
md" # Tentative merging it all"

# ╔═╡ d0b197de-8b32-47f8-b95f-56d3979015fd
Tevap_b_out = substitute(T2B, x=>0.0)

# ╔═╡ 44cb764c-ec69-4396-8fbb-6889847d4d09
Tbh_in = T1_in

# ╔═╡ 938b5761-2a8e-430f-a280-f0eae4e68bad
eq_in_bh = Tbh_in ~ Tevap_b_out

# ╔═╡ b7a8f43b-1232-4e55-8b5e-e7eab13f933c
Tbh_out = substitute(sT2, x=>0.0)

# ╔═╡ 2d779ec5-2be7-4079-b8c0-e16ccab558fb
Tevap_b_in = T2iB

# ╔═╡ 87b6739c-9f07-4962-8123-fe6cbf3d9de4
eq_out_bh = Tbh_out ~ Tevap_b_in

# ╔═╡ bb715483-8993-44d7-b793-c2ce72ee55b6
eq_flow1 = m1 ~ m2

# ╔═╡ 4199640a-375d-4c96-936a-14cddc2b92ed
eq_flow2 = ṁ1 ~ m2

# ╔═╡ 4ff535c4-444e-4893-ac72-83925eb70c4f
@named overall = NonlinearSystem([eq_cont_2, eq_sat_T1, eq_cont_1, eq_length, eq_in_bh, eq_out_bh, eq_flow1, eq_flow2], [LA, LB, ṁ2, m2, m1, T2iA, T2iB, T1_in], [L, Tsat, U_A, ṁ1, U_B, h1iA, hdew, T2iB, T1oB, cp1, cp2, P, Ts, UA, Ue1, Ue2, Lbh])

# ╔═╡ c578bed8-02a7-4a05-9765-00142a85aa9e
ov = complete(simplify(overall))

# ╔═╡ a9587932-391a-46ad-84e0-9c4760e534af
md"## Utilities"

# ╔═╡ bfc9c576-5670-4a51-a811-b073e260d205
T_range_evap = -20:0.1:20.0

# ╔═╡ 6a121d29-491f-4b7b-8295-367e6683add8
md"""Borehole wall temperature: **Tw = $(@bind Twall_ Scrubbable(T_range_evap, default = 10.0)) °C**"""

# ╔═╡ d273ecb6-6fb6-4fdb-a820-faec866ace7a
md"""Fluid freezing temperature: **Tf = $(@bind Tfreeze_ Scrubbable(T_range_evap)) °C**"""

# ╔═╡ f5873723-097c-49d2-90de-eba8d34b3b27
md"""
Evaporation temperature: $(@bind Tevap_ Scrubbable(T_range_evap; default=-5.0)) °C
"""

# ╔═╡ 9f09862a-7bcc-47fb-89f9-d3ab4748e0f0
χ_range_ = fluid != "Water" ? (0.0:1.0:100.0) : NaN

# ╔═╡ f15060e9-b543-464b-b97d-287afe15efc4
md"""Concentration: $(@bind concentration Scrubbable(χ_range_)) % (mass)"""

# ╔═╡ 1ba60683-8370-40c9-abf0-d4a58713d708
Pevap_ = PropsSI("P", "T", 273.15+Tevap_, "Q", 1.0, refrigerant)

# ╔═╡ 984f3e2a-c965-45fb-a492-c64de4d3e3d3
md"""
Evaporation pressure: **$(round(Pevap_/1e5, digits = 2))** bar
"""

# ╔═╡ c550217f-e4be-4a4a-931e-39d100294441
h_dew_val = PropsSI("HMASS", "T", 273.15+Tevap_, "Q", 1.0, refrigerant)

# ╔═╡ 4e2df2e4-da65-483b-be8f-479e60048f30
md"""
Dew enthalpy evaporator: **$(round(h_dew_val/1e3))** kJ/kg
"""

# ╔═╡ e5862511-38f4-40a5-baa2-7b9a440b398c
h_comp_val = PropsSI("HMASS", "T", 273.15+Tevap_+ΔT_SH_, "P", Pevap_, refrigerant)

# ╔═╡ 9f562f5f-6fb9-4e45-a019-1c85f9910189
md"""
Enthalpy suction compressor: **$(round(h_comp_val/1e3))** kJ/kg
"""

# ╔═╡ 3ec0a698-d433-4022-a42e-30158f951f5a
cp_dew_val = PropsSI("C", "T", 273.15+Tevap_, "Q", 1.0, refrigerant)

# ╔═╡ eadefc6c-7ba5-464d-99bc-d9055134a721
md"""
Dew thermal capacity evaporator: **$(round(cp_dew_val))** J/kg/K
"""

# ╔═╡ a0bf0d50-1125-4c34-b60e-4505e9edee95
T_range_cond = 25:0.1:65.0

# ╔═╡ ea5d77e6-c42d-436c-9941-6c54e6d45ec8
md"""
Condensation temperature: $(@bind Tcond_ Scrubbable(T_range_cond; default=40.0)) °C
"""

# ╔═╡ 7ab03084-d007-4657-a5cd-74b13c671011
Pcond_ = PropsSI("P", "T", 273.15+Tcond_, "Q", 1.0, refrigerant)

# ╔═╡ 08db51d7-89ee-4a27-a83f-627b5fbd108a
md"""
Condensation pressure: **$(round(Pcond_/1e5, digits = 2))** bar
"""

# ╔═╡ c9885140-68b4-4ab9-a56e-0934bc54de6e
Tbrine_ = Tevap_+ΔT_SH_+ΔT_brine_

# ╔═╡ c9635ea5-8aa6-447e-9590-f6a90e053e6f
md"""
Inlet brine temperature: **$(Tbrine_) °C**
"""

# ╔═╡ ff754da1-688e-40d7-99d1-a91b6a8b4ba2
if fluid != "Water"
	fluid_call = "INCOMP::"*fluid*"[$(concentration/100)]"
else
	fluid_call = fluid
end

# ╔═╡ 0f7c8016-d127-4842-ac74-273f3a5ccc79
cp_brine_ = Unitful.ustrip(PropsSI("C", "T", Tbrine_+273.15, "Q", 0.0, fluid_call))

# ╔═╡ a7eb03ee-b55f-41fc-82b3-7f81a7ec9db0
md"""
Fluid thermal capacity: cp = **$(round(cp_brine_)) J/kg/K**
"""

# ╔═╡ caf0077d-73ff-417f-a711-f8c1479c167d
ps = Dict(
	L => L_,
	Tsat => Tevap_,
	U_A => U_A_,
	U_B => U_B_,
	h1iA => 350e3, 
	hdew => h_dew_val,
	T1oB => Tevap_ + ΔT_SH_,
	T2iB => Tbrine_,
	cp1 => cp_dew_val,
	cp2 => cp_brine_,
	P => P_,

	LA => L/2.0,
	LB => L/2.0,
	ṁ1 => ṁ_r_,
	ṁ2 => 1.1,
	T2iA => Tsat,
)

# ╔═╡ 00b587d0-d6cc-4513-8ef7-b12ac5d08b13
prob_tz = NonlinearProblem(tz, [], ps, symbolic_u0 = true)

# ╔═╡ 390f1bf1-5083-463b-bda3-8efb5e300aba
sol = NonlinearSolve.solve(prob_tz)

# ╔═╡ 440558b1-ef48-47c3-9ca2-b13616be666c
md"""
Fluid mass flow: **$(round(sol[ṁ2], digits = 2)) kg/s**
"""

# ╔═╡ 16530d22-5691-4b89-b3d3-9975cf67aad1
sol[LA]

# ╔═╡ 072008d4-f6ea-4b3e-8594-f8e15da7c576
sol[LB]

# ╔═╡ 425a6b31-f45a-4a8c-8dfa-c862c37c7470
sol[ṁ2]

# ╔═╡ 16097dcb-4d8e-47ea-abbf-0bc6f3a8e1e7
sol[T2iA]

# ╔═╡ a140827a-d88e-4f6d-82e5-9f1902912e3e
sol.ps[T1oB]

# ╔═╡ 337fd039-79ec-4159-86b8-ce77424492f8
sol[T2iA]

# ╔═╡ 63e66cb4-688a-4bc1-b9f2-a82866188471
T1A_(x) = ps[Tsat]

# ╔═╡ e7d6c2fb-a428-40a5-9b5e-dbfd356975da
T2A_(x) = T2A_n(x, sol[LA], sol[T2iA], ps[Tsat], ps[P], ps[cp2], ps[U_A], sol[ṁ2])

# ╔═╡ ab8b13eb-fb11-4945-a925-6f04e3f2c6cf
T2A_(1.0)

# ╔═╡ d2c870bd-f8a5-4f13-a457-6ed981b6d4af
T1B_(x) = T1B_n(x, sol[LB], ps[T1oB], ps[T2iB], ps[P], ps[cp1], ps[cp2], ps[U_B], ps[ṁ1], sol[ṁ2])

# ╔═╡ 3b0515ea-f6b2-4501-b081-50a8131056e4
T1B_(1.0)

# ╔═╡ 7152dfbf-d16f-4255-a03e-1e5e88e80c9f
T1B_(sol[LB])

# ╔═╡ 5ec817be-3af6-49bc-b993-36cfa1e815da
T2B_(x) = T2B_n(x, sol[LB], ps[T1oB], ps[T2iB], ps[P], ps[cp1], ps[cp2], ps[U_B], ps[ṁ1], sol[ṁ2])

# ╔═╡ 75dc9775-fc40-4986-af2c-3b29b9678d8c
T2B_(sol[LB])

# ╔═╡ d56115df-acd5-44de-9a38-affd71731ea5
begin
	T1_n=build_function(sT1, xbh, Lbh, UA, Ue1, Ue2, m1, m2, cp1, cp2, T1_in, Ts, expression = Val{false})
	T1_(x) =T1_n(x,Lbh_,Upp_,Upw_,Upw_,sol[ṁ2], sol[ṁ2], cp_brine_, cp_brine_, T2A_(0.0), Twall_)
end

# ╔═╡ 9f474167-459f-4bb0-a566-47254a587852
begin
	T2_n=build_function(sT2, xbh, Lbh, UA, Ue1, Ue2, m1, m2, cp1, cp2, T1_in, Ts, expression = Val{false})
	T2_(x) =T2_n(x,Lbh_,Upp_,Upw_,Upw_,sol[ṁ2], sol[ṁ2], cp_brine_, cp_brine_, T2A_(0.0), Twall_)
end

# ╔═╡ a27c5fc7-dad2-44f3-ab65-5736fefc874d
begin
    # Create a figure with two subplots side by side
    fig = Figure(size = (1200, 600))
	y_range = (-10, 15)
    
    # First subplot: (e.g., the "bh" figure)
    ax1 = Axis(fig[1, 1], xlabel="x (m)", ylabel="Temperature (°C)", title="Borehole",xreversed=true)
    ylims!(ax1, y_range)
    lines!(ax1, 0..Lbh_, x -> T1_(x), color=:blue)
    lines!(ax1, 0..Lbh_, x -> T2_(x), color=:blue)
    hlines!(ax1, Twall_, color=:yellow)
    
    # Second subplot: (e.g., the counterflow HX)
    ax2 = Axis(fig[1, 2], xlabel="x (m)", ylabel="Temperature (°C)", title="Evaporator")
    ylims!(ax2, y_range)
    lines!(ax2, sol[LA]..ps[L], x -> T1B_(x - sol[LA]), color=:green, linewidth=2, label="Refrigerant")
    lines!(ax2, sol[LA]..ps[L], x -> T2B_(x - sol[LA]), color=:blue, linewidth=2, label="Brine")
    lines!(ax2, 0..sol[LA], x -> T1A_(x), color=:green, linewidth=2, label="Refrigerant")
    lines!(ax2, 0..sol[LA], x -> T2A_(x), color=:blue, linewidth=2, label="Brine")
	hlines!(ax2, T2B_(ps[L] - sol[LA]), color=:blue, linestyle =:dash)
    
    fig
end


# ╔═╡ 1be8ef60-8f28-44fd-9309-a89d4eb2c0e3
ps_ov = Dict(
	L => L_,
	Lbh => Lbh_,
	Tsat => Tevap_,
	U_A => U_A_,
	U_B => U_B_,
	h1iA => 350e3, 
	hdew => h_dew_val,
	T1oB => Tevap_ + ΔT_SH_,
	T2iB => Tbrine_,
	cp1 => cp_dew_val,
	cp2 => cp_brine_,
	P => P_,
	Ts => Twall_,
	UA => U_A_,
	Ue1 => U_A_,
	Ue2 => U_A_,

	LA => L/2.0,
	LB => L/2.0,
	ṁ1 => ṁ_r_,
	ṁ2 => 1.1,
	m1 => 1.1,
	m2 => 1.1,
	T2iA => Tsat,
	T2iB => Tsat,
	T1_in => Twall_,
)

# ╔═╡ f099df79-afed-40dc-8996-688188d34ac7
prob_ov = NonlinearProblem(ov, [], ps_ov, symbolic_u0 = true)

# ╔═╡ 247dad59-682d-4f99-997b-fb3d67cf0fdd
sol_ov = NonlinearSolve.solve(prob_ov)

# ╔═╡ 276bcba0-9a4d-442f-9d78-acb64bb284d9
λ_brine_ = Unitful.ustrip(PropsSI("L", "T", Tbrine_+273.15, "Q", 0.0, fluid_call))

# ╔═╡ f4ca4219-2e3a-412a-9c06-98b72cb63ae0
md"""
Fluid thermal conductivity: λ = **$(round(λ_brine_, digits=2)) W/m/K**
"""

# ╔═╡ 7ed12177-a1e4-45c7-aa0e-bd28fa633479
ρ_brine_ = Unitful.ustrip(PropsSI("D", "T", Tbrine_+273.15, "Q", 0.0, fluid_call))

# ╔═╡ bea8fea1-7248-4fcc-892c-6f8f3c6e46cb
md"""
Fluid density: **ρ = $(round(ρ_brine_, digits=1)) kg/m³**
"""

# ╔═╡ f853312b-5da5-43d3-9b51-51dec3bce19f
μ_brine_ = Unitful.ustrip(PropsSI("V", "T", Tbrine_+273.15, "Q", 0.0, fluid_call))

# ╔═╡ 5184eed4-0e01-4d05-aefb-84b6b7811ddb
md"""
Fluid dynamic viscosity: μ = **$(round(μ_brine_, digits=3)) Pa⋅s**
"""

# ╔═╡ 6514d469-f5b0-4567-b808-433066cca548
md"---"

# ╔═╡ c70fe3b1-ad7d-4538-8487-dd7c0c1296a2
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(30px, 5%);
    	padding-right: max(500px, 5%);
	}
</style>
"""

# ╔═╡ Cell order:
# ╟─8b583195-b96d-4cd6-8bb5-f09512c38ea6
# ╟─a27c5fc7-dad2-44f3-ab65-5736fefc874d
# ╟─4b6d717a-edb6-4292-bcb7-5dc44e560637
# ╟─7e0fdfd5-dc29-4cb1-9ee5-82442ec8933c
# ╟─81a7c7ba-9cba-4b50-98c8-e14f5edf64d3
# ╟─ac5f77da-a6b7-4be1-8fb6-523c10f2a147
# ╟─9e09fbba-ecdf-4615-8baf-f25a7cafb195
# ╟─b754c0a6-4579-40f2-b03c-31142ef5bb2f
# ╟─2dffa4af-cb4e-42f0-be0e-417734d119d4
# ╟─6a121d29-491f-4b7b-8295-367e6683add8
# ╟─0a284427-2f1e-4d6f-b3d9-32b0b54c4564
# ╟─7182187c-c0a3-4b8b-bd7e-06718192666a
# ╟─d273ecb6-6fb6-4fdb-a820-faec866ace7a
# ╟─f15060e9-b543-464b-b97d-287afe15efc4
# ╟─a7eb03ee-b55f-41fc-82b3-7f81a7ec9db0
# ╟─f4ca4219-2e3a-412a-9c06-98b72cb63ae0
# ╟─5184eed4-0e01-4d05-aefb-84b6b7811ddb
# ╟─bea8fea1-7248-4fcc-892c-6f8f3c6e46cb
# ╟─4af116ad-3eac-479d-8c44-2e90ab2021c2
# ╟─c9635ea5-8aa6-447e-9590-f6a90e053e6f
# ╟─89c842e7-73ac-4258-83e4-a98ca23e5a12
# ╟─937426c1-8075-45a7-9d92-8a729b1ac251
# ╟─436b93eb-fc29-4529-9e6d-2aa16e808f1c
# ╟─64756cb4-776f-4512-8f2f-fc11b9715859
# ╟─8e6ffc58-6d02-4707-8504-4c917281e29f
# ╟─30d089ba-e823-408c-9073-50b3f37a1c4c
# ╟─4937015c-8461-4cfe-8c83-51c3544f3abf
# ╟─440558b1-ef48-47c3-9ca2-b13616be666c
# ╟─25b4a88c-5c58-4d67-8d0f-ecbad5f5afae
# ╟─3fc0ed88-9d4d-45c3-a0f2-48e38e9caa43
# ╟─3703ecf4-8a7d-44f7-a59a-f429ad768205
# ╟─f5873723-097c-49d2-90de-eba8d34b3b27
# ╟─984f3e2a-c965-45fb-a492-c64de4d3e3d3
# ╟─30e8753a-b57c-4382-910d-fcaf08bbe6b3
# ╟─4e2df2e4-da65-483b-be8f-479e60048f30
# ╟─eadefc6c-7ba5-464d-99bc-d9055134a721
# ╟─9f562f5f-6fb9-4e45-a019-1c85f9910189
# ╟─4bf19bc9-29f9-493f-830c-a6beb30fc8d6
# ╟─ea5d77e6-c42d-436c-9941-6c54e6d45ec8
# ╟─08db51d7-89ee-4a27-a83f-627b5fbd108a
# ╟─ed6c8661-293b-4a95-a714-1fcf498883fa
# ╟─1b08d63a-db57-44c7-b1ec-10acbcaa7951
# ╠═d4aba512-dd8c-49a5-a7b1-eccf7ca24eab
# ╠═9902671a-a681-4b14-9e07-424eb673eef1
# ╠═d5cbb87a-d35d-4aff-b213-c89ed193619e
# ╠═c7bd0055-8fbe-4861-9284-1e08d89bcad3
# ╠═1f23f14b-69e4-4366-acd8-8ded3a6ca847
# ╠═b6d31712-03e2-4981-95c4-d3621f9d0ed6
# ╠═7dd35a9e-12d1-4728-bae4-49c909cd7902
# ╠═8781d78b-1308-496f-8e80-aa1cd334dd12
# ╟─39ca3785-1a2f-4f1f-8877-76c9841ce098
# ╠═7274b23c-f22e-4b2a-86b1-5976d25c363c
# ╠═b146988c-40cb-4576-97f1-32950538af16
# ╠═ca92e844-5285-47de-a83d-5ae5c872e44b
# ╠═2a4fca86-f5f9-4e0d-8a89-8898ea190cca
# ╠═af00af29-d180-425b-88d5-92a9c8d86b79
# ╠═08d860a7-c75a-4fed-a752-63894211aee6
# ╠═be4aa932-54b7-4017-9c06-59104b11aea0
# ╠═cc07e1f9-8877-43e8-99c3-2d293bc00427
# ╟─1481a416-efa8-4ef9-9866-c159ce0f3e75
# ╠═92665c1f-eb3a-4c13-9f8a-5eed60e54835
# ╠═ee30d6f9-b902-4d87-9ce7-f245573642ee
# ╟─c87f7981-ccf7-4b22-ada3-7e7e944f1bf4
# ╠═744924ad-810b-4bd0-9417-165f6a48ed9a
# ╠═e89934b8-49a3-4e32-8bb1-fdc6021163fa
# ╠═007afd7d-7941-4bda-99a4-fe02d6b1d0be
# ╠═3f1d4587-27fe-458a-8765-7465d29d6688
# ╟─470dccc9-6eca-4275-a532-d03889791b06
# ╠═79558ba6-4469-4a40-966b-490f6d8e20f0
# ╠═78d75ba8-0c50-48ff-8c4c-1f75498eb109
# ╠═14ee0cbd-81e0-4da0-8b81-8020d8c0175e
# ╠═5b8e04d3-0080-4755-9f48-526871e0e434
# ╠═84adc973-3aad-4ed1-aab3-4ff3142da7d8
# ╠═e8bcc3d7-fc8f-4960-b63b-90e769502f0e
# ╠═caf0077d-73ff-417f-a711-f8c1479c167d
# ╠═00b587d0-d6cc-4513-8ef7-b12ac5d08b13
# ╠═390f1bf1-5083-463b-bda3-8efb5e300aba
# ╠═16530d22-5691-4b89-b3d3-9975cf67aad1
# ╠═072008d4-f6ea-4b3e-8594-f8e15da7c576
# ╠═425a6b31-f45a-4a8c-8dfa-c862c37c7470
# ╠═16097dcb-4d8e-47ea-abbf-0bc6f3a8e1e7
# ╠═a140827a-d88e-4f6d-82e5-9f1902912e3e
# ╠═63e66cb4-688a-4bc1-b9f2-a82866188471
# ╠═4656a6dc-e314-4110-8d8c-e20ece6417bd
# ╠═b0f357cc-f0ad-4a28-912f-9ee0517d49bc
# ╠═337fd039-79ec-4159-86b8-ce77424492f8
# ╠═e7d6c2fb-a428-40a5-9b5e-dbfd356975da
# ╠═ab8b13eb-fb11-4945-a925-6f04e3f2c6cf
# ╠═84ed2576-78fc-4779-a964-85bfbde1e01c
# ╠═39120525-093c-494e-8494-da3bd900b650
# ╠═d2c870bd-f8a5-4f13-a457-6ed981b6d4af
# ╠═3b0515ea-f6b2-4501-b081-50a8131056e4
# ╠═ba97eec7-b2fe-4c9c-8b8d-c6f20cde17b3
# ╠═4fe5b696-dc7f-4d52-b4e7-afbe9e9d7ca6
# ╠═5ec817be-3af6-49bc-b993-36cfa1e815da
# ╠═75dc9775-fc40-4986-af2c-3b29b9678d8c
# ╠═7152dfbf-d16f-4255-a03e-1e5e88e80c9f
# ╟─d77ec1b5-ebbb-463d-ae71-9cd880e5dd1a
# ╠═c6acaf2a-5f42-4932-acb1-632ad13689e5
# ╠═cea2d602-810a-4dc4-a0fe-0fb76d680b58
# ╠═0610a6d3-7b33-415c-bf55-9c16750785e7
# ╠═e8887d73-60e3-430a-b998-e864e7e4bcf3
# ╠═5bf5942e-46b4-4716-8e57-062a54d436e0
# ╠═3f01b035-8fb7-4031-8a11-ef745c83606b
# ╠═56e7aae4-89ae-4bb9-90a4-1ebbd6300d12
# ╠═97fd0f42-4e1e-453f-a239-925003eb103c
# ╠═1c648223-2e60-4d35-83a8-d5488f79884b
# ╠═12b47c63-28d1-4a78-a733-23e10dbfe214
# ╠═690e2521-eda6-4037-a7c3-69c48f82a25a
# ╠═6b85b101-74e2-4d27-9c27-1fff81802c53
# ╠═15752913-4b37-4de4-9337-86ace6dcf64f
# ╠═3d776a5a-a56f-4094-a15e-e2e275f8ffe5
# ╠═72c892e4-64fe-44b5-bcfb-0673289be93f
# ╠═20b6adea-53c7-4d6c-91d4-9784fd6d8113
# ╠═ef63e77a-f931-4704-9774-dee7f37ee9eb
# ╠═56016b33-a3c6-4f3b-a6e9-37687349d984
# ╠═39f5e272-fe0b-456b-a3f3-d268d1f2756a
# ╠═40a43ccc-1fad-440a-b033-66061b4dffa9
# ╠═1bc06242-a9e0-4836-a061-4aa896eb11b7
# ╠═8d53b90c-124f-433a-a9d3-87c83fddfa22
# ╠═2b9a45fb-05b1-4b54-9c96-88289d79aee8
# ╠═01ad16df-5b68-40a0-9e5e-16fea88ebc8e
# ╠═77caced1-75ee-42ba-8483-59cf3d10fa14
# ╠═043640e4-5b45-4e91-a0e1-c286be1b17c1
# ╠═1f43adb4-7812-45e6-b4fb-a89f207cceb9
# ╠═840b2e1c-922a-4664-a622-05a4da45a3ca
# ╠═3b34cf1d-3390-4ad5-9e69-b56dca416dbd
# ╟─d56115df-acd5-44de-9a38-affd71731ea5
# ╟─9f474167-459f-4bb0-a566-47254a587852
# ╟─b7f99ca5-5531-4bfb-ab80-0b8513578d63
# ╠═d0b197de-8b32-47f8-b95f-56d3979015fd
# ╠═44cb764c-ec69-4396-8fbb-6889847d4d09
# ╠═938b5761-2a8e-430f-a280-f0eae4e68bad
# ╠═b7a8f43b-1232-4e55-8b5e-e7eab13f933c
# ╠═2d779ec5-2be7-4079-b8c0-e16ccab558fb
# ╠═87b6739c-9f07-4962-8123-fe6cbf3d9de4
# ╠═bb715483-8993-44d7-b793-c2ce72ee55b6
# ╠═4199640a-375d-4c96-936a-14cddc2b92ed
# ╠═4ff535c4-444e-4893-ac72-83925eb70c4f
# ╠═c578bed8-02a7-4a05-9765-00142a85aa9e
# ╠═1be8ef60-8f28-44fd-9309-a89d4eb2c0e3
# ╠═f099df79-afed-40dc-8996-688188d34ac7
# ╠═247dad59-682d-4f99-997b-fb3d67cf0fdd
# ╟─a9587932-391a-46ad-84e0-9c4760e534af
# ╟─bfc9c576-5670-4a51-a811-b073e260d205
# ╠═9f09862a-7bcc-47fb-89f9-d3ab4748e0f0
# ╠═1ba60683-8370-40c9-abf0-d4a58713d708
# ╠═7ab03084-d007-4657-a5cd-74b13c671011
# ╠═c550217f-e4be-4a4a-931e-39d100294441
# ╠═e5862511-38f4-40a5-baa2-7b9a440b398c
# ╠═3ec0a698-d433-4022-a42e-30158f951f5a
# ╠═a0bf0d50-1125-4c34-b60e-4505e9edee95
# ╠═c9885140-68b4-4ab9-a56e-0934bc54de6e
# ╠═ff754da1-688e-40d7-99d1-a91b6a8b4ba2
# ╠═0f7c8016-d127-4842-ac74-273f3a5ccc79
# ╠═276bcba0-9a4d-442f-9d78-acb64bb284d9
# ╠═7ed12177-a1e4-45c7-aa0e-bd28fa633479
# ╠═f853312b-5da5-43d3-9b51-51dec3bce19f
# ╟─6514d469-f5b0-4567-b808-433066cca548
# ╠═bc6d9b90-734e-4d9d-be0f-104a3b99aea2
# ╠═8857ee4a-9620-4071-ae78-6e63ea0b7f0e
# ╠═a9eddcc7-ac94-483e-8622-e6c9dca8351b
# ╠═0405ca05-0ca9-4b6c-ac93-91e5f50cab91
# ╠═c691da35-ff83-4dc6-af9f-2be03d95187a
# ╠═f6b9843f-bb5e-4275-bea2-05e9d594207f
# ╠═1c61fab6-f56b-460b-870f-8e3a83337312
# ╠═fbf3c417-8b31-4256-8981-84b0f0ce34b0
# ╠═fa00c042-0b93-4635-ac75-4513e61671ec
# ╠═c70fe3b1-ad7d-4538-8487-dd7c0c1296a2
