"""
example: uses Finite State Projection to solve the burr08 model
"""

from cmepy.models import burr08

import numpy
import cmepy.recorder
import cmepy.fsp.solver
import cmepy.fsp.simple_expander
import cmepy.domain
import cmepy.state_enum
import cmepy.statistics

import fsp_example_util

def main():
    """
    Solve burr08 model using experimental dense FSP

    This only uses the simple expansion routine ...
    """

    model = burr08.create_model()

    # construct the initial domain states, and a StateEnum
    # defining a bijection between them and array indices
    initial_states = cmepy.domain.from_iter((model.initial_state, ))
    initial_state_enum = cmepy.state_enum.StateEnum(initial_states)
    # define the initial dense probability with respect to this state enum
    dense_p_0 = initial_state_enum.pack_distribution({model.initial_state : 1.0})

    expander = cmepy.fsp.simple_expander.SimpleExpander(
        model.transitions,
        depth = 3,
    )

    fsp_solver = cmepy.fsp.solver.create(
        model,
        initial_states,
        initial_state_enum, # nb state enum required!
        expander,
        time_dependencies = burr08.create_time_dependencies(),
        p_0 = dense_p_0,
    )

    # define time steps:
    # this problem is initially stiff so
    # we begin with some finer time steps
    # before changing to coarser steps
    time_steps = numpy.concatenate((
        numpy.linspace(0.0, 1.0, 10),
        numpy.linspace(2.0, 16.0, 15)
    ))

    # we want the error of the solution at the
    # final time to be bounded by epsilon
    epsilon = 1.0e-2
    num_steps = numpy.size(time_steps)
    # define how much error is tolerated per step
    max_error_per_step = epsilon / num_steps

    # create recorder to record species counts
    recorder = cmepy.recorder.create(
        (model.species, model.species_counts)
    )

    domains = []

    for i, t in enumerate(time_steps):
        print 'STEP t = %g' % t
        fsp_solver.step(t, max_error_per_step)
        if i % 3 == 0:
            # Get the dense distribution as an array.
            # This is defined with respect to the
            # current bijection between array indices
            # and states, which is represented as a
            # cmepy.state_enum.StateEnum instance. You
            # can find that in fsp_solver.solver.solver.domain_enum
            p_dense, p_sink = fsp_solver.y

            # The Recorder expects a sparse distribution,
            # so we have to do the sparse conversion explicitly.
            domain_enum = fsp_solver.solver.solver.domain_enum
            p_sparse = domain_enum.unpack_distribution(p_dense)

            # record sparse distribution
            print 'recording solution and domain'
            recorder.write(t, p_sparse)

            # also store a copy of the domain so we can plot it later
            domains.append(numpy.array(fsp_solver.domain_states))
    print 'OK'

    print 'plotting solution and domain'
    fsp_example_util.plot_solution_and_domain(
        recorder[('A', 'B')],
        domains
    )

if __name__ == '__main__':
    main()

