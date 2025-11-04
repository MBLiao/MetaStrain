"""Provide minimization of metabolic adjustment (MOMA)."""

from typing import TYPE_CHECKING, Optional

from optlang.symbolics import Zero, add

from cobra.util import solver as sutil
import pandas as pd
from cobra.flux_analysis.parsimonious import pfba


if TYPE_CHECKING:
    from cobra.core import Model, Solution


def moma(model, reference_fluxes, linear: bool = False) -> "Solution":

    with model:
        add_moma(model=model, reference_fluxes=reference_fluxes, linear=linear)
        solution = model.optimize()
    return solution


def add_moma(model, reference_fluxes, linear: bool = False) -> None:

    if "moma_old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for MOMA.")

    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)

    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="moma_old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)

    obj_vars = []
    for protID in reference_fluxes.index:
        r = model.reactions.get_by_id(protID)
        flux = reference_fluxes[protID]

        if linear:
            components = sutil.add_absolute_expression(
                model,
                r.flux_expression,
                name="moma_dist_" + r.id,
                difference=flux,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(
                r.flux_expression - dist,
                lb=flux,
                ub=flux,
                name="moma_constraint_" + r.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist**2)

    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)




