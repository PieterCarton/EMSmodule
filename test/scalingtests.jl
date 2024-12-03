# Examples from https://jump.dev/JuMP.jl/stable/tutorials/getting_started/tolerances/
using JuMP
import HiGHS
import SCS
n, ε = 13, 0.0234
N = 2^n
model = Model(SCS.Optimizer)
@variable(model, x[1:N] >= 0)
@objective(model, Min, x[1])
@constraint(model, sum(x) == 1)
z = [(-1)^((i & (1 << j)) >> j) for j in 0:n-1, i in 0:N-1]
@constraint(model, z * x .>= 1 - ε)
set_attribute(model, "eps_abs", 1e-5)
set_attribute(model, "eps_rel", 1e-5)

optimize!(model)
is_solved_and_feasible(model)
value(x[1])
1 - n * ε / 2
report = primal_feasibility_report(model)
maximum(values(report))
violated_variables = filter(xi -> value(xi) < 0, x)
y = first(violated_variables)
value(y)
report[LowerBoundRef(y)]